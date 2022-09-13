# Alignement d'embedding par programmation dynamique

## Setup your environment

Clone the repository:

```bash
git clone 
```

Move to the new directory:

```bash
cd alignment_embedding
```

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Install [mamba](https://github.com/mamba-org/mamba):

```bash
conda install mamba -n base -c conda-forge
```

Create the `embedding` conda environment:
```
mamba env create -f binder/environment.yml
```

Load the `embedding` conda environment:
```
conda activate embedding
```

To deactivate an active environment, use

```
conda deactivate
```
## Modify the configuration file

`[files]` section
put file name of prot to align (emb + fasta)
`[align]` section 
choose True or False for each alignment

## Run code

```
python path/to/src/main.py
```