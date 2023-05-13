# Metagenomes metageneration metatool

>  A simple and smart tool for generating Illumina readings of metagenomic communities with a given properties

### Installation

To get the tool clone the git repository:

```bash
https://github.com/nvaulin/Metagenomes_Generations_Tool.git
cd Metagenomes_Generations_Tool
```

Create a `conda/mamba` environment with necessary packages and activate it:

```bash
conda env create -f environment.yml
conda activate metageneration
```

### Usage


To run the script, just call it from the directory where the tools is located:

```bash
python Metagenome_generation.py [PHENOTYPE] ...
```

To perform the test run use the `2_species` phenotype:
```bash
python Metagenome_generation.py 2_species
```

With the real baseline phenotypes its better to select `n` core species with the `ncore` option.
To test the pathways correction use the `example_pathways.txt` file:
```bash
python Metagenome_generation.py -p Health -c 10 --pathways example_pathways.txt
```

To get more information about the particular script, run:

```bash
python Metagenome_generagtion.py  --help
```


### Uninstallation

To uninstall the tool remove the conda environment and delete the cloned folder:
```python
conda remove --name metageneration --all
rm -rf Metagenomes_Generations_Tool
```

### Citation

If you use these tool, please cite as:
- Vaulin N., Chechenina A., Ivanov. A, Uliantsev V. (2023). Metagenomes metageneration metatool. Природа, 2013, 443, 27–28
