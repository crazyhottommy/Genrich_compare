
# ATACseq

### download data for comparing

I downloaded the same data sets used in this paper https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz533/5519166

six SRA accession numbers SRR891269-SRR891274. There are three biological replicates generated using 50 000 cells per replicate and other three generated using 500 cells per replicate.

type them in https://ewels.github.io/sra-explorer/

```bash
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891273/SRR891273_1.fastq.gz -o SRR891273_GSM1155962_GM12878_ATACseq_500_Rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891273/SRR891273_2.fastq.gz -o SRR891273_GSM1155962_GM12878_ATACseq_500_Rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891271/SRR891271_1.fastq.gz -o SRR891271_GSM1155960_GM12878_ATACseq_50k_Rep4_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891271/SRR891271_2.fastq.gz -o SRR891271_GSM1155960_GM12878_ATACseq_50k_Rep4_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891270/SRR891270_1.fastq.gz -o SRR891270_GSM1155959_GM12878_ATACseq_50k_Rep3_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891270/SRR891270_2.fastq.gz -o SRR891270_GSM1155959_GM12878_ATACseq_50k_Rep3_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_1.fastq.gz -o SRR891269_GSM1155958_GM12878_ATACseq_50k_Rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_2.fastq.gz -o SRR891269_GSM1155958_GM12878_ATACseq_50k_Rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891274/SRR891274_1.fastq.gz -o SRR891274_GSM1155963_GM12878_ATACseq_500_Rep3_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891274/SRR891274_2.fastq.gz -o SRR891274_GSM1155963_GM12878_ATACseq_500_Rep3_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891272/SRR891272_1.fastq.gz -o SRR891272_GSM1155961_GM12878_ATACseq_500_Rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891272/SRR891272_2.fastq.gz -o SRR891272_GSM1155961_GM12878_ATACseq_500_Rep1_Homo_sapiens_OTHER_2.fastq.gz
```


```bash

git clone 

#make a meta file with 4 columns
cat meta.txt
sample  fastq factor    read
SRR891269       SRR891269_GSM1155958_GM12878_ATACseq_50k_Rep2_Homo_sapiens_OTHER_1.fastq.gz     atac    R1
SRR891269       SRR891269_GSM1155958_GM12878_ATACseq_50k_Rep2_Homo_sapiens_OTHER_2.fastq.gz     atac    R2
SRR891270       SRR891270_GSM1155959_GM12878_ATACseq_50k_Rep3_Homo_sapiens_OTHER_1.fastq.gz     atac    R1
SRR891270       SRR891270_GSM1155959_GM12878_ATACseq_50k_Rep3_Homo_sapiens_OTHER_2.fastq.gz     atac    R2
SRR891271       SRR891271_GSM1155960_GM12878_ATACseq_50k_Rep4_Homo_sapiens_OTHER_1.fastq.gz     atac    R1
SRR891271       SRR891271_GSM1155960_GM12878_ATACseq_50k_Rep4_Homo_sapiens_OTHER_2.fastq.gz     atac    R2
SRR891272       SRR891272_GSM1155961_GM12878_ATACseq_500_Rep1_Homo_sapiens_OTHER_1.fastq.gz     atac    R1
SRR891272       SRR891272_GSM1155961_GM12878_ATACseq_500_Rep1_Homo_sapiens_OTHER_2.fastq.gz     atac    R2
SRR891273       SRR891273_GSM1155962_GM12878_ATACseq_500_Rep2_Homo_sapiens_OTHER_1.fastq.gz     atac    R1
SRR891273       SRR891273_GSM1155962_GM12878_ATACseq_500_Rep2_Homo_sapiens_OTHER_2.fastq.gz     atac    R2
SRR891274       SRR891274_GSM1155963_GM12878_ATACseq_500_Rep3_Homo_sapiens_OTHER_1.fastq.gz     atac    R1
SRR891274       SRR891274_GSM1155963_GM12878_ATACseq_500_Rep3_Homo_sapiens_OTHER_2.fastq.gz     atac    R2

python sample2json.py /n/scratchlfs/informatics/mtang/Genrich_compare/atac/ /n/scratchlfs/informatics/mtang/Genrich_compare/atac/meta.txt

# a samples.json file will be generated
cat samples.json | head 
{
    "SRR891269": {
        "atac": {
            "R1": [
                "/n/scratchlfs/informatics/mtang/Genrich_compare/atac/SRR891269_GSM1155958_GM12878_ATACseq_50k_Rep2_Homo_sapiens_OTHER_1.fastq.gz"
            ],
            "R2": [
                "/n/scratchlfs/informatics/mtang/Genrich_compare/atac/SRR891269_GSM1155958_GM12878_ATACseq_50k_Rep2_Homo_sapiens_OTHER_2.fastq.gz"
            ]
        }


nano config.yaml
```
this is ATACseq data, there is no control
set  "control: False"


### download reference genome 

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

## get a bed file for streches of Ns in the genome and provide it to Genrich -E

python ~/apps/Genrich/findNs.py hg38.fa hg38_Ns.bed

## make a bwa index
bwa index hg39.fa 
```

