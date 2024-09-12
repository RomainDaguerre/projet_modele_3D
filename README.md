# Evaluation Quality 3D (EQ3D)

PROGRAMME D’EVALUATION DE LA QUALITE D’UN MODELE 3D DE
PROTEINE
Programme qui évalue la qualité d'un modèle 3D de protéine à partir des potentiels statistiques DOPE (https://www.dsimb.inserm.fr/~gelly/data/dope.par)

<p align="center">
    <img alt="Made with Python" src="https://img.shields.io/badge/Made%20with-Python-1f425f.svg?color=%23539fc9">
</p>

## Installation

### Clone the project:

#### SSH:
```bash
git clone git@github.com:RomainDaguerre/projet_court.git
```
OR
#### HTTPS:
```bash
git clone https://github.com/RomainDaguerre/projet_court
```

### Move to the new directory:
```bash
cd projet_court
```

### Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Create a Conda environment

```bash
conda env create -f environment.yml
```

### Activate the Conda environment

```bash
conda activate environment
```

### You can also update the conda environment with:

```bash
mamba env update -f binder/environment.yml
```

### To deactivate an conda active environment, use

```
conda deactivate
```

## Usage as command line tool

### Run EQ3D on a test file:
```bash
python EQ3D.py test/1bta.pdb
```

### Run EQ3D:
```bash
python EQ3D.py your_pdb_file.pdb
```

### Run EQ3D_group:
For run a group of pdb file.
```bash
python EQ3D_group.py your_file_where_pdb/
```


