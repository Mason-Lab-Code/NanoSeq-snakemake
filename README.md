# NanoSeq-snakemake
Description...

### Install

### Set up initial directory structure

### Link FASTQ files with compatible file name system

### Download reference files

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
