  # BLAST-SNAKE 

### BLAST-Snake will take any number of query sequences in FASTA format and
###  output a plotted phylogeny.

#### To install dependencies to run BLAST-Snake, perform the following steps:

#### Install conda.
  #### 1. Download a conda version from https://docs.conda.io/en/latest/miniconda.html
  #### 2. Run the downloaded file from the /home/ directory.
             $ bash xxxxx-latest-Linux-x86_64.sh
  #### 3. Accept defaults and install.

#### Create environment.
  #### 1. Navigate to the BLAST-Snake directory.
  #### 2. Create the environment from file.
             $ conda env create --file envs/blast-snake.yaml
  #### 3. Activate the new environment
             $ conda activate blast-snake


#### To run BLAST-Snake:

  #### Navigate to the BLAST-snake/ directory.
  #### Ensure environment is activated.
  #### Place query FASTAs in data/query/
  #### Ensure queries have the suffix ".fasta"

  #### To run the Snake for all files in data/query/
				     $ snakemake

  #### To run the Snake for a specific file in data/query/
				     $ snakemake <specific file>
