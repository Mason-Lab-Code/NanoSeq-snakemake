# NanoSeq-snakemake

This is a Snakemake workflow for the complete NanoSeq bioinformatics pipeline (https://github.com/cancerit/NanoSeq), designed to run on the University of York's HPC cluster, Viking. 

### Install



### Set up initial directory structure

Move to the directory in which you would like to execute the workflow. Create subdirectories 00_raw/ and 01_ref/. 
```
mkdir 00_raw 01_ref
```
### Link FASTQ files with compatible file name system

FASTQ file names must have the following structure: SAMPLE_CONDITION_read1.fastq.gz and SAMPLE_CONDITION_read1.fastq.gz (e.g. donor1_treated_read1.fastq.gz)
You don't need to change the names of the original files. Instead, create symlinks to the original files, using the new naming structure, inside the 00_raw/ subdirectory. 
```
ln -s /path/to/original-read1.fastq.gz 00_raw/SAMPLE_CONDITION_read1.fastq.gz
ln -s /path/to/original-read2.fastq.gz 00_raw/SAMPLE_CONDITION_read2.fastq.gz
```

### Download reference files

Download the required reference files into the 01_ref/ subdirectory. 

```

```

### Create config.yaml

### Set up tmux session

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

```
snakemake --slurm --default-resources slurm_account=bio-cancerinf-2020 slurm_partition=nodes --jobs 24 --use-conda --conda-frontend conda
```
