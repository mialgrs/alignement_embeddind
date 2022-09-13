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
## Execute the alignement script

### Move to the right directory

Use the command `cd` to change directory.
```bash
cd src
```
### Modify the configuration file
 - `[files]` section

Put the name of the fasta and t5emb file of the 2 protein.

 - `[align]` section 

Choose which alignment you need by indicating `True` or `False` for each.

### Run the script

```
python main.py
```