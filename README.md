# Genrich_compare
snakemake workflow comparing Genrich and MACS2

### workflow of the pipeline


### Dependencies

* [bwa/0.7.15](https://github.com/lh3/bwa) for aligning the reads.

* [Genrich](https://github.com/jsh58/Genrich)

git clone the repo on July 22th

```bash
git clone https://github.com/jsh58/Genrich
cd Genrich
make

```

* [MACS2](https://github.com/taoliu/MACS)

macs2 2.1.2

```bash
conda install macs2
```
### How to run it

```bash

ssh bio1

## start a screen session
screen

git clone https://github.com/crazyhottommy/Genrich_compare

## make samples.json file
python sample2json.py /path/to/fastq/dir meta.txt

```

The `config.yaml` file contains various parameters one need to change.
Now it can handle single-end or pair-end fastqs; ChIPseq with or without Input/IgG controls
and ATACseq data.

```bash
## whether it is from fastq files or from bam files
from_fastq: True

## the reads are paired end or single end
paired_end: True


# >70bp long, <70bp false
long_reads: True

# what's the control for IP
# if it is ChIP-seq data
# this is the same as the third colmn of the meta.txt factor column

control: 'Input'

## if it is ChIP-seq data without Input/IgG control or ATACseq data (without control by nature)
## set control : False


## the reference fasta path
ref_fa: "/n/regal/informatics_public/reference_genome_by_tommy/ucsc/mm9/mm9.fa"

macs2_g: mm
macs2_pvalue: 1e-5

# other macs2 arguments
macs2_args: -f BAM

Genrich_args: "-r -p 0.05 -a 200 "

#number of reads downsample to, I set to 50 million, if reads number smaller than
## 50 million, downsample will keep the orignal reads
target_reads: 50000000
```


```bash
## set up the config.yaml file

nano config.yaml


## dry run
snakemake -np

## real run
snakemake -j 24

## submit to slurm
./pyflow_Genrich_compare.sh
```

### To do

- [ ] Adatpor trimming
- [ ] Add benchmark directive for logging time and memory usage.
- [ ] Genrich with replicates https://github.com/jsh58/Genrich#multiple-replicates
- [ ] IDR framework https://github.com/kundajelab/idr compare it with Genrich. IDR is fast.
- [ ] conda env and docker container 
