# Answer to skills_test
A snakemake pipeline for annotation and filtering VCF files.

## List of used software tools:
* python =3.8.5
* R =4.0.3
* java =11.0.9.1
* picard =2.23.3
* vcftools =0.1.16
* bcftools =1.10.2
* plink2 =2.00a2.3LM
* tabix =1.10.2-3
* VariantAnnotation r package =1.36.0
* snakemake =5.10.0

## Operating system requirements
* moreutils =0.63-1

## Compilation instructions
* No compilation needed.

## How to run the workflow

### Method 1

#### Step 1: Configure workflow
Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 2: Execute workflow
Test your configuration by performing a dry-run via

    $ snakemake --use-conda -n


#### Step 3: Execute the workflow locally via

    $ snakemake --use-conda --cores $N

### Method 2

#### Step 1: Configure workflow
Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 2: Creating an environment using environment.yaml

    $ conda env create --name smake_env --file envs/environment.yaml

#### Step 3: Activating the environment

	$ conda activate smake_env

#### Step 4: Running workflow

	$ snakemake

To exit the environment, execute

    $ conda deactivate
