import os
import argparse
import yaml
import subprocess
import pandas as pd
from Bio import Entrez
from pysat.examples.hitman import Hitman

from scripts.Genomes_db_update import update_genomes, write_multifasta

GENOMES_DIR = 'genomes'
RESULTS_DIR = 'results'


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
    parser.add_argument('-t', '--threads', default=1, help='number of threads (cores)')
    parser.add_argument('email', default='example@email.com', nargs='?',
                        help='Email address for Entrez requests')
    parser.add_argument('api_key', default=None, nargs='?',
                        help='NCBI API key for Entrez requests (if any)')

    return parser.parse_args()


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
    print(pathways_input)
    if os.path.isfile(pathways_input):
        with open(pathways_input, 'r') as f:
            pathways = f.readlines()
    elif ',' in pathways_input:
        pathways = pathways_input.split(',')
    else:
        raise ValueError('Invalid input. Please provide a path to a file or a comma-separated string.')
    return pathways


if __name__ == '__main__':
    pheno = parse_args().phenotype
    metagenome_file = parse_args().metagenome_file
    pathways = parse_args().pathways
    n_core = parse_args().n_core
    n_threads = parse_args().threads
    email = parse_args().email
    api_key = parse_args().api_key

    Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key

    os.makedirs(RESULTS_DIR, exist_ok=True)

    abundances = pd.read_csv(os.path.join('baseline_phenotypes', pheno + '.tsv'), sep='\t', header=None)
    abundances.rename({0: 'species', 1: 'abundance'}, axis=1, inplace=True)
    if n_core:
        n_core = min(n_core, len(abundances))
        abundances = abundances.sort_values(by='abundance').head(n_core)
    pathways_db = pd.read_csv(os.path.join('Databases',
                                           'Pathways_MetaCyc.txt'), sep='\t').dropna(inplace=True)
    # pathways_db = filter_pathways_db(pathways_db)
    # pathways_db = preprocess_pathways_db(pathways_db)
    if pathways is not None:
        print('Reading required pathways...')
        pathways_specified = read_pathways(pathways)
        species_to_refill = find_minimal_refill(abundances[0].to_list(),
                                                pathways_specified)
        abundances = append_species_refill(abundances, species_to_refill)
    prepared_abudances = update_genomes(GENOMES_DIR, abundances)
    wr_code = write_multifasta(prepared_abudances, GENOMES_DIR)
    print('\n')

    iss_params = {
        '-g': os.path.join(GENOMES_DIR, 'multifasta.fna'),
        '--abundance_file': os.path.join(RESULTS_DIR, 'abundances.txt'),
        '-m': 'miseq',
        '-o': os.path.join(RESULTS_DIR, 'miseq_reads'),
        '--cpus': n_threads
    }
    with open('iss_params.yml', 'r') as f:
        yaml_params = yaml.safe_load(f)
        iss_params = iss_params | yaml_params

    iss_cmd = ['iss', 'generate'] + [str(item) for pair in iss_params.items() for item in pair]
    result = subprocess.run(iss_cmd)
    if result.returncode == 0:
        print('\nThe metagenome was successfully generated!')
    else:
        print('\nThe metagenome generation completed with errors.')
