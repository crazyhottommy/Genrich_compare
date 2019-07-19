# Genrich_compare
snakemake workflow comparing Genrich and MACS2

### workflow of the pipeline


### Dependencies

### How to run it

```bash

ssh bio1

## start a screen session
screen

git clone https://github.com/crazyhottommy/Genrich_compare

## make samples.json file
python sample2json.py /path/to/fastq/dir meta.txt

## dry run
snakemake -np

## real run
snakemake -j 24

## submit to slurm
./pyflow_Genrich_compare.sh
```
