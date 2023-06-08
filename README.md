# NanoSeq-snakemake

This is a Snakemake workflow for the complete NanoSeq bioinformatics pipeline (https://github.com/cancerit/NanoSeq), designed to run on the University of York's HPC cluster, Viking. 

### Install

To use this Snakemake workflow, simply move to the directory in which you would like to execute the workflow, and clone this repository.  
```
git clone https://github.com/Mason-Lab-Code/NanoSeq-snakemake/
```

### Set up initial directory structure

Create subdirectories 00_raw/ and 01_ref/. 
```
mkdir 00_raw 01_ref
```
### Link FASTQ files with compatible file name system

FASTQ file names must have the following structure: SAMPLE_CONDITION_read1.fastq.gz and SAMPLE_CONDITION_read1.fastq.gz (e.g. donor1_treatment1_read1.fastq.gz)
You don't need to change the names of the original files. Instead, create symlinks to the original files, using the new naming structure, inside the 00_raw/ subdirectory. 
```
ln -s </path/to/original-read1.fastq.gz> <00_raw/SAMPLE_CONDITION_read1.fastq.gz>
ln -s </path/to/original-read2.fastq.gz> <00_raw/SAMPLE_CONDITION_read2.fastq.gz>
```

### Download reference files

Download the required reference files into the 01_ref/ subdirectory. 

```
cd 01_ref
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
```
You will also need genomic mask BED files provided in this Google folder: https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL
Genomic mask files required: SNP.sorted.GRCh38.bed.gz, SNP.sorted.GRCh38.bed.gz.tbi, NOISE.sorted.GRCh38.bed.gz, NOISE.sorted.GRCh38.bed.gz.tbi
Download these locally and move them to the 01_ref/ subdirectory on Viking (throuhg WinSCP or similar). 

### Create config.yaml

You will need to edit the config.yaml file to specify the names of your files. 

Under "SAMPLES:" enter the names of your samples, which should be the first part of your FASTQ file names (before the first underscore). 

Under "TYPES:" enter all the conditions of your experiment, which should be the second part of your FASTQ files names (between the first and second underscores). 

Under "DUPLEX_TYPES:" enter the conditions of your experiment that are duplex. 

Under "UNDILUTED_TYPES:" enter the conditions of your experiment that are undiluted (matched normals). 

The config.yaml file should look something like this. Make sure you indent as below, and include the colon at the end of each line. 

```
SAMPLES:
    donor1:
    donor2:

TYPES:
    treatment1:
    treatment2:

DUPLEX_TYPES:
    treatment1:
    treatment2:

UNDILUTED_TYPES:
    control:
```

### Set up tmux session

This is an overview on how to use tmux, which is a terminal multiplexer that allows you to invoke Snakemake, and allow it to continue running even when you disconnect from Viking. 

```
module load tools/tmux
tmux new -s session_name
```
Type Ctrl+B then D to detach from the session (it will still be active) and go back to the main terminal. 

You can now exit the terminal and log out of Viking, and the tmux session will still be active and accessible when you log back in. 

List out active tmux sessions - you should see the session that you just created. 
```
tmux ls
```
Attach to the tmux session 
```
tmux attach-session -t session_name
```
Kill the tmux session after the snakemake workflow has completed. 
```
tmux kill-session -t session_name
```
### Dry run 

Invoke a "dry run" of the workflow to ensure everything is set up properly. You should see all of the tasks that Snakemake 
```
module load tools/snakemake
snakemake -n
```
Create a directed acyclic graph to visualise the steps that snakemake will run. 
```
snakemake --dag | dot -Tsvg > dag.svg
```

### Run snakemake 

Start a new tmux session (as above) and invoke Snakemake to start executing the workflow for real. Snakemake should start submitting jobs to Viking. 

```
# Create tmux session
# Attach to tmux session
# Run Snakemake
snakemake --slurm --default-resources slurm_account=bio-cancerinf-2020 slurm_partition=nodes --jobs 24 --use-conda --conda-frontend conda
# Detach from tmux session
```
Check on the progress of submitted jobs using squeue. 
```
squeue -u <username>
```

### Making sense of outputs

Path to output VCF files (SNVs): 10_analysis/<SAMPLE_CONDITION-vs-UNDILUTED>/tmpNanoSeq/post/results.muts.vcf.gz
Path to output VCF files (indels): 10_analysis/<SAMPLE_CONDITION-vs-UNDILUTED>/tmpNanoSeq/post/results.indel.vcf.gz
Path to useful summary of variants called: 10_analysis/<SAMPLE_CONDITION-vs-UNDILUTED>/tmpNanoSeq/post/summary.txt
Path to mutational profile visualisations: 10_analysis/<SAMPLE_CONDITION-vs-UNDILUTED>/results.trinuc-profiles.pdf
Path to contamination check output: 11_contamination_check/<SAMPLE_CONDITION.out> (explanation here: https://github.com/cancerit/NanoSeq)
Path to efficiency estimation output: 12_efficiency_estiamte/<SAMPLE_CONDITION.tsv> (explanation here: https://github.com/cancerit/NanoSeq)
