import os
import urllib
import argparse
import gzip
import urllib.request
import pandas as pd
from Bio import Entrez, SeqIO, Seq
from pysat.examples.hitman import Hitman


def parse_args():
    parser = argparse.ArgumentParser(
        usage='Metagenome_generation.py [PHENO]  ...',
        description='''Generate Illumina reads for the  described metagenome in the FILE. '''
    )

    parser.add_argument('-p', '--phenotype', default='2_species', nargs='?',
                        help='the base phenotype for metagenome construction ("Health", "HIV")')
    parser.add_argument('-m', '--metagenome_file', default=None, nargs='?',
                        help='read metagenome composition from the file (tsv with species and abudances)')
    parser.add_argument('pathways', default=None, nargs='?',
                        help='read matebolic pathways to account from the file (each pathway on the new line')
    parser.add_argument('n_core', default=None, nargs='?',
                        help='number of core species to leave in metagenome')
    parser.add_argument('email', default='example@email.com', nargs='?',
                        help='Email address for Entrez requests')
    parser.add_argument('api_key', default=None, nargs='?',
                        help='NCBI API key for Entrez requests (if any)')

    return parser.parse_args()


def generate_ncbi_search_terms(specie):
    terms_for_search = []
    specie_for_search = f'"{specie}"' + '[Organism]'
    filter_for_relevant = '(all[filter] NOT anomalous[filter])'
    filter_for_latest = 'latest[filter]'
    filter_for_refseq = '"latest refseq"[filter]'
    filter_for_genbank = '"latest genbank"[filter]'
    filters_for_taxonomy = ['" AND taxonomy check ok"[filter]', '']
    filters_for_complete = ['" AND complete genome"[filter]', '']
    filters_for_db = [filter_for_refseq, filter_for_genbank, filter_for_latest]

    for filter_for_complete in filters_for_complete:
        for filter_for_taxonomy in filters_for_taxonomy:
            for filter_for_db in filters_for_db:
                terms_temp = [specie_for_search, filter_for_relevant, filter_for_db]
                terms_for_search.append(' AND '.join(terms_temp) + filter_for_taxonomy + filter_for_complete)
    return specie_for_search, terms_for_search


def update_genomes(genomes_dir, baseline_abudances):
    if not os.path.isdir(genomes_dir):
        os.makedirs(genomes_dir)
    prepared_abudances = pd.DataFrame(columns=['taxid', 'abudance'])
    assembly_summary_cols = ['RsUid', 'GbUid', 'AssemblyAccession', 'LastMajorReleaseAccession',
                             'LatestAccession', 'ChainId', 'AssemblyName', 'UCSCName', 'EnsemblName',
                             'Taxid', 'Organism', 'SpeciesTaxid', 'SpeciesName', 'AssemblyType',
                             'AssemblyStatus', 'AssemblyStatusSort', 'WGS', 'GB_BioProjects',
                             'GB_Projects', 'RS_BioProjects', 'RS_Projects', 'BioSampleAccn',
                             'BioSampleId', 'Biosource', 'Coverage', 'PartialGenomeRepresentation',
                             'Primary', 'AssemblyDescription', 'ReleaseLevel', 'ReleaseType',
                             'AsmReleaseDate_GenBank', 'AsmReleaseDate_RefSeq', 'SeqReleaseDate',
                             'AsmUpdateDate', 'SubmissionDate', 'LastUpdateDate',
                             'SubmitterOrganization', 'RefSeq_category', 'AnomalousList',
                             'ExclFromRefSeq', 'PropertyList', 'FromType', 'Synonym', 'ContigN50',
                             'ScaffoldN50', 'AnnotRptUrl', 'FtpPath_GenBank', 'FtpPath_RefSeq',
                             'FtpPath_Assembly_rpt', 'FtpPath_Stats_rpt', 'FtpPath_Regions_rpt',
                             'Busco', 'SortOrder', 'Meta']

    assembly_status_translation = {'Complete Genome': 1,
                                   'Scaffold': 2,
                                   'Contig': 3}

    for entry in range(len(baseline_abudances)):
        specie, abudance = baseline_abudances.iloc[entry]
        assemblies_summary = pd.DataFrame(columns=assembly_summary_cols)
        specie_for_search, terms_for_search = generate_ncbi_search_terms(specie)
        tax_id = Entrez.read(Entrez.esearch(db="taxonomy", term=specie))['IdList'][0]
        fna_filename = str(tax_id) + '.fna.gz'
        if fna_filename not in os.listdir(genomes_dir):
            for term in terms_for_search:
                assembly_ids = Entrez.read(Entrez.esearch(db="assembly", term=term, retmax=100))['IdList']
                if len(assembly_ids) > 0:
                    break
            for assembly_id in assembly_ids:
                assembly_summary = Entrez.read(Entrez.esummary(db="assembly", id=assembly_id, report="full"))
                assembly_summary = pd.DataFrame(assembly_summary['DocumentSummarySet']['DocumentSummary'],
                                                columns=assembly_summary_cols)
                assemblies_summary = pd.concat([assemblies_summary, assembly_summary])

            assemblies_summary.replace({"AssemblyStatus": assembly_status_translation}, inplace=True)
            assemblies_sorted = assemblies_summary.sort_values(['AssemblyStatus', 'ScaffoldN50',
                                                                'Coverage', 'LastUpdateDate', 'ContigN50'],
                                                               ascending=[True] + [False] * 4)
            links = assemblies_sorted[['FtpPath_RefSeq', 'FtpPath_GenBank']].head(1)
            if links.FtpPath_RefSeq.values[0]:
                url = links.FtpPath_RefSeq.values[0]
            else:
                url = links.FtpPath_GenBank.values[0]
            label = os.path.basename(url)
            url = url.removeprefix('ftp://')
            url = 'https://' + url + '/' + label + '_genomic.fna.gz'
            urllib.request.urlretrieve(url, os.path.join(genomes_dir, fna_filename))

            specie_abudance_df = pd.DataFrame([tax_id, abudance])
            prepared_abudances = pd.concat([prepared_abudances, specie_abudance_df])

    prepared_abudances.abudance = prepared_abudances.abudance / prepared_abudances.abudance.sum()
    return prepared_abudances


def write_multifasta(prepared_abudances, genomes_dir):
    sequences = []
    for tax_id in prepared_abudances.taxid:
        seq_record, seq_str = Seq(), ''
        seq_record.id = tax_id
        seq_record.name = tax_id
        with gzip.open(os.path.join(genomes_dir, f'{tax_id}.fna.gz'), "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_str += record
        seq_record.seq = seq_str
        sequences.append(seq_record)

    wr_code = SeqIO.write(sequences, os.path.join(genomes_dir, 'multifasta.fna'), "fasta")
    return wr_code


def filter_pathways_db(pathways_db):
    junk_tax_ranges = ['Metazoa', 'Embryophyta', 'Tracheophyta', 'Pinidae', 'Brassicales', 'Gunneridae',
                       'Spermatophyta', 'Vertebrata <vertebrates>', 'Fungi', 'Eukaryota', 'Mammalia',
                       'cellular organisms', 'Viridiplantae', 'Magnoliopsida', 'Fungi // Viridiplantae',
                       'Fungi // Metazoa', 'Viruses // Metazoa']
    for junk_entry in junk_tax_ranges:
        pathways_db = pathways_db[pathways_db['Taxonomic-Range'] != junk_entry]
    selection = pathways_db['Taxonomic-Range'].str.contains('bact')
    return pathways_db[selection]


def preprocess_pathways_db(pathways_db):
    pathways_db.Pathways = pathways_db[['Pathways']].replace('\(?<i>.+</i>\)?', '', regex=True)
    pathways_db['Species'] = pathways_db['Species'].str.split('//').apply(
        lambda l: [x.strip() for x in l]).values.tolist()
    all_species = [item.strip() for sublist in pathways_db['Species'].values.tolist() for item in sublist]
    pathways_db_exp = pathways_db[['Pathways', 'Species']].explode('Species')
    euks = ['Homo sapiens', 'Arabidopsis thaliana', 'Saccharomyces cerevisiae', 'Glycine max', 'Pisum sativum',
            'Rattus norvegicus', 'Solanum lycopersicum', 'Oryza sativa', 'Nicotiana tabacum']
    mask = pathways_db_exp.Species.isin(euks)
    pathways_db_exp = pathways_db_exp[~mask]
    pathways_db_ct = pd.crosstab(index=pathways_db_exp.Pathways, columns=pathways_db_exp.Species)
    pathways_db_ct = pathways_db_ct.astype(bool).astype(int)
    return pathways_db_ct


def do_hits(metagenome: list, metabolic_needs: list[list]):
    h = Hitman(solver='m22', htype='lbx')
    metagenome = set(metagenome)
    needs_to_hit = []
    for metabolic_need in metabolic_needs:
        metabolic_need = set(metabolic_need)
        if not set.intersection(metabolic_need, metagenome):
            needs_to_hit.append(metabolic_need)
    for need in needs_to_hit:
        h.hit(need)
    return h.get()


def find_minimal_refill(metagenome, metabolites_specified, pathways_db):
    cols = pathways_db.columns
    selected_bathways = pathways_db.loc[metabolites_specified].astype('bool')
    metabolic_needs = selected_bathways.apply(lambda x: list(cols[x.values]), axis=1).to_list()
    return do_hits(metagenome, metabolic_needs)


def append_species_refill(abudances, species_to_refill):
    abundance_level = abudances[1].mean()
    abundances_refill = pd.DataFrame([species_to_refill,
                                      [abundance_level] * len(species_to_refill)],
                                     index=[0, 1]).transpose()
    abudances_new = pd.concat([abudances, abundances_refill])
    abudances_new[1] = abudances_new[1] / abudances_new[1].sum()
    return abudances_new

def read_pathways(pathways_input):
    if os.path.isfile(pathways_input.input):
        with open(pathways_input, 'r') as f:
            pathways = f.readlines()
    elif ',' in pathways_input.input:
        pathways = pathways_input.input.split(',')
    else:
        raise ValueError('Invalid input. Please provide a path to a file or a comma-separated string.')
    return pathways

if __name__ == '__main__':
    pheno = parse_args().phenotype
    metagenome_file = parse_args().metagenome_file
    pathways = parse_args().pathways
    n_core = parse_args().n_core
    email = parse_args().email
    api_key = parse_args().api_key

    Entrez.email = email
    Entrez.api_key = api_key

    genomes_dir = os.path.join('genomes')
    baseline_abudances = pd.read_csv(os.path.join('baseline_phenotypes', pheno + '.tsv'), sep='\t', header=None)
    pathways_db = pd.read_csv(os.path.join('Metabolites_database',
                                           'Pathways_MetaCyc.txt'), sep='\t').dropna(inplace=True)
    pathways_db = filter_pathways_db(pathways_db)
    pathways_db = preprocess_pathways_db(pathways_db)

    pathways_specified = read_pathways(pathways)
    species_to_refill = find_minimal_refill(baseline_abudances[0].to_list(),
                                            pathways_specified)
    abudances_new = append_species_refill(baseline_abudances, species_to_refill)
    prepared_abudances = update_genomes(genomes_dir, abudances_new)
    wr_code = write_multifasta(prepared_abudances, genomes_dir)
