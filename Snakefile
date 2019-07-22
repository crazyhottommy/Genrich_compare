shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

localrules: all


CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os

SAMPLES = sorted(FILES.keys())

## list all samples by sample_name and sample_type
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        MARK_SAMPLES.append(sample + "_" + sample_type)

ALL_SAMPLES = list(set(MARK_SAMPLES))

# which sample_type is used as control for calling peaks: e.g. Input, IgG...
# or set control: False if there is no control
CONTROL = config["control"]
if CONTROL:
    CONTROLS = [sample for sample in MARK_SAMPLES if CONTROL in sample]
    CASES = [sample for sample in MARK_SAMPLES if CONTROL not in sample]
else:
    CASES = MARK_SAMPLES

ALL_BAM = expand("03aln/{sample}.sorted.bam", sample = ALL_SAMPLES)

## peaks and bigwigs
ALL_MACS_PEAKS = []
ALL_GENRICH_PEAKS = []

for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    if CONTROL:
        control = sample + "_" + CONTROL
        if control in CONTROLS:
            ALL_MACS_PEAKS.append("09peak_macs2/{}_vs_{}_macs2_peaks.xls".format(case, control))
            ALL_GENRICH_PEAKS.append("09peak_Genrich/{}_vs_{}_Genrich_peaks.xls".format(case, control))
    else:
        ALL_MACS_PEAKS.append("09peak_macs2/{}_macs2_peaks.xls".format(case))
        ALL_GENRICH_PEAKS.append("09peak_Genrich/{}_Genrich_peaks.bed".format(case))


ALL_DOWNSAMPLE_BAM = expand("04aln_downsample/{sample}-downsample.sorted.bam", sample = ALL_SAMPLES)
ALL_INDEX = expand("03aln/{sample}.sorted.bam.bai", sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand("04aln_downsample/{sample}-downsample.sorted.bam.bai", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("03aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_BIGWIG = expand("07bigwig/{sample}.bw", sample = ALL_SAMPLES)
ALL_QC = ["10multiQC/multiQC_log.html"]


TARGETS = []
TARGETS.extend(ALL_BAM)
TARGETS.extend(ALL_DOWNSAMPLE_BAM)
TARGETS.extend(ALL_INDEX)
TARGETS.extend(ALL_MACS_PEAKS)
TARGETS.extend(ALL_GENRICH_PEAKS)

rule all:
    input: TARGETS

## get a list of fastq.gz files for the same mark, same sample
def get_input_files(wildcards):
    sample_name = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    if not config["paired_end"]:
        return FILES[sample_name][mark]['R1']
    ## start with bam files
    if not config["from_fastq"]:
        return FILES[sample_name][mark]

def get_R1_files(wildcards):
    sample_name = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    if config["from_fastq"]:
        if config["paired_end"]:
            return  FILES[sample_name][mark]['R1']

def get_R2_files(wildcards):
    sample_name = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    if config["from_fastq"]:
        if config["paired_end"]:
            return  FILES[sample_name][mark]['R2']

## when there are multiple fastq.gz for the same sample, same mark, merge them.
if config["paired_end"]:
    rule merge_fastqs:
        input:
            r1 = get_R1_files,
            r2 = get_R2_files
        output:
            "01seq/{sample}_R1.fastq.gz", "01seq/{sample}_R2.fastq.gz"
        log: "00log/{sample}_merge_fastq.log"
        params:
            jobname = "{sample}"
        threads: 1
        message: "merging fastqs {input}: {threads} threads"
        shell:
            """
            gunzip -c {input.r1} | gzip > {output[0]} 2> {log}
            gunzip -c {input.r2} | gzip > {output[1]} 2>> {log}
            """
else:
    rule merge_fastqs:
        input: get_input_files
        output: temp("01seq/{sample}.fastq.gz")
        log: "00log/{sample}_unzip"
        threads: CLUSTER["merge_fastqs"]["n"]
        params: jobname = "{sample}"
        message: "merging fastqs gunzip -c {input} > {output}"
        shell: "gunzip -c {input} | gzip > {output} 2> {log}"


        # get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and dupliated reads by samblaster -r
# samblaster should run before samtools sort

## align from fastq files
## deal with paired end and single end data
## see https://bitbucket.org/snakemake/snakemake/pull-requests/148/examples-for-paired-read-data-and/diff
## some hints:
## https://bitbucket.org/snakemake/snakemake/issues/37/add-complex-conditional-file-dependency
## conditional rules?
## https://groups.google.com/forum/#!msg/snakemake/qX7RfXDTDe4/cKZBfc_PAAAJ


## when feeding @RG to bwa mem,  literal tab will cause the resulting bam file header violating the
## SAM specification at @PG line. https://github.com/lh3/bwa/issues/83. This will interfere with IGV and picard
## to process the bam files.
## and https://github.com/Duke-GCB/bespin-cwl/issues/6
## BWA after version BWA-0.7.16a (r1181) with verbose >=1 will  this problem.
## now feed bwa with read group by escaping the tab
## Before: -R "@RG\tID:1\tLB:LIBRARY\tPL:illumina\tSM:sample1\tPU:AB1234"
## Now: -R "@RG\\tID:1\\tLB:LIBRARY\\tPL:illumina\\tSM:sample1\\tPU:AB1234"


if config["from_fastq"] and config["paired_end"]:
    rule align:
        input:  r1 = "01seq/{sample}_R1.fastq.gz",
                r2 = "01seq/{sample}_R2.fastq.gz"
        output: "03aln/{sample}.sorted.bam", "03aln/{sample}.sorted.bam.bai", "00log/{sample}.align"
        threads: CLUSTER["align"]["n"]
        params:
                jobname = "{sample}",
                ## add read group for bwa mem mapping, change accordingly if you know PL:ILLUMINA, LB:library1 PI:200 etc...
                rg = "@RG\\tID:{sample}\\tSM:{sample}"
        message: "aligning bwa {input}: {threads} threads"
        log:
            bwa = "00log/{sample}.align",
            markdup = "00log/{sample}.markdup"
        run:
            if config["long_reads"]:
                shell(
                    r"""
                    bwa mem -t 5 -M -v 1 -R '{params.rg}' {config[ref_fa]} {input[0]} {input[1]} 2> {log.bwa} \
                    | samblaster 2> {log.markdup} \
                    | samtools view -Sb -F 4 - \
                    | samtools sort -m 2G -@ 5 -T {output[0]}.tmp -o {output[0]}
                    samtools index {output[0]}
                    """)
            ## short reads < 70bp
            ## Probably one of the most important is how many mismatches you will allow between a read and a potential mapping location for that location to be considered a match.
            ## The default is 4% of the read length, but you can set this to be either another proportion of the read length, or a fixed integer
            else:
                shell(
                    r"""
                    bwa aln -t 5 {config[ref_fa]} {input[0]} 2> {log.bwa} > 03aln/{wildcards.sample}_R1.sai
                    bwa aln -t 5 {config[ref_fa]} {input[1]} 2>> {log.bwa} > 03aln/{wildcards.sample}_R2.sai
                    bwa sampe -r '{params.rg}' {config[ref_fa]} 01aln/{wildcards.sample}_R1.sai 01aln/{wildcards.sample}_R2.sai {input[0]} {input[1]} 2>> {log.bwa} \
                    | samblaster 2> {log.markdup} \
                    | samtools view -Sb -F 4 - \
                    | samtools sort -m 2G -@ 5 -T {output[0]}.tmp -o {output[0]}
                    rm 03aln/{wildcards.sample}_R1.sai 03aln/{wildcards.sample}_R2.sai
                    samtools index {output[0]}
                    """)


if config["from_fastq"] and not config["paired_end"]:
    rule align:
        input:  "01seq/{sample}.fastq.gz"
        output: "03aln/{sample}.sorted.bam", "03aln/{sample}.sorted.bam.bai", "00log/{sample}.align"
        threads: CLUSTER["align"]["n"]
        params:
                jobname = "{sample}",
                ## add read group for bwa mem mapping, change accordingly if you know PL:ILLUMINA, LB:library1 PI:200 etc...
                rg = "@RG\\tID:{sample}\\tSM:{sample}"
        message: "aligning bwa {input}: {threads} threads"
        log:
            bwa = "00log/{sample}.align",
            markdup = "00log/{sample}.markdup"
        run:
            if config["long_reads"]:
                shell(
                    r"""
                    bwa mem -t 5 -M -v 1 -R '{params.rg}' {config[ref_fa]} {input} 2> {log.bwa} \
                    | samblaster -r 2> {log.markdup} \
                    | samtools view -Sb -F 4 - \
                    | samtools sort -m 2G -@ 5 -T {output[0]}.tmp -o {output[0]}
                    samtools index {output[0]}
                    """)
            ## short reads < 70bp
            else:
                shell(
                    r"""
                    bwa aln -t 5 {config[ref_fa]} {input} 2> {log.bwa} > 01aln/{wildcards.sample}.sai
                    bwa samse -r '{params.rg}' {config[ref_fa]} 01aln/{wildcards.sample}.sai {input} 2>> {log.bwa} \
                    | samblaster -r 2> {log.markdup} \
                    | samtools view -Sb -F 4 - \
                    | samtools sort -m 2G -@ 5 -T {output[0]}.tmp -o {output[0]}
                    rm 01aln/{wildcards.sample}.sai
                    samtools index {output[0]}
                    """)

rule flagstat_bam:
    input:  "03aln/{sample}.sorted.bam"
    output: "03aln/{sample}.sorted.bam.flagstat"
    log:    "00log/{sample}.flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule down_sample:
    input: "03aln/{sample}.sorted.bam", "03aln/{sample}.sorted.bam.bai", "03aln/{sample}.sorted.bam.flagstat"
    output: "04aln_downsample/{sample}-downsample.sorted.bam", "04aln_downsample/{sample}-downsample.sorted.bam.bai"
    log: "00log/{sample}_downsample.log"
    threads: 5
    params: jobname = "{sample}"
    message: "downsampling for {input}"
    run:
        import re
        import subprocess
        with open (input[2], "r") as f:
            # fifth line contains the number of mapped reads
            line = f.readlines()[4]
            match_number = re.match(r'(\d.+) \+.+', line)
            total_reads = int(match_number.group(1))

        target_reads = config["target_reads"] # 15million reads  by default, set up in the config.yaml file
        if total_reads > target_reads:
            down_rate = target_reads/total_reads
        else:
            down_rate = 1

        shell("sambamba view -f bam -t 5 --subsampling-seed=3 -s {rate} {inbam} | samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} 2> {log}".format(rate = down_rate, inbam = input[0], outbam = output[0], log = log))

        shell("samtools index {outbam}".format(outbam = output[0]))

rule make_bigwigs:
    input : "04aln_downsample/{sample}-downsample.sorted.bam", "04aln_downsample/{sample}-downsample.sorted.bam.bai"
    output: "07bigwig/{sample}.bw"
    log: "00log/{sample}.makebw"
    threads: 5
    params: jobname = "{sample}"
    message: "making bigwig for {input}"
    shell:
        """
        bamCoverage -b {input[0]} --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 5 --extendReads 200 -o {output} 2> {log}
        """

## there is an Input or IgG control
if CONTROL:
    rule call_peaks_macs2:
        input: control = "04aln_downsample/{control}-downsample.sorted.bam", case="04aln_downsample/{case}-downsample.sorted.bam"
        output: bed = "09peak_macs2/{case}_vs_{control}_macs2_peaks.xls"
        log: "00log/{case}_vs_{control}_call_peaks_macs2.log"
        params:
            name = "{case}_vs_{control}_macs2",
            jobname = "{case}_vs{control}_macs2",
            custom = config.get("macs2_args")
        message: "call_peaks macs2 {input}: {threads} threads"
        shell:
            """
           ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
            source activate macs2
            macs2 callpeak -t {input.case} \
                -c {input.control} -g {config[macs2_g]} \
                --outdir 09peak_macs2 -n {params.name} -p {config[macs2_pvalue]} {params.custom} &> {log}
            """
else:
    rule call_peaks_macs2:
        input: case = "04aln_downsample/{case}-downsample.sorted.bam"
        output: bed = "09peak_macs2/{case}_macs2_peaks.xls"
        log: "00log/{case}_call_peaks_macs2.log"
        params:
            name = "{case}_macs2",
            jobname = "{case}_macs2",
            custom = config.get("macs2_args")
        message: "call_peaks macs2 {input}: {threads} threads"
        shell:
            """
           ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
            source activate macs2
            macs2 callpeak -t {input.case} \
                -g {config[macs2_g]} \
                --outdir 09peak_macs2 -n {params.name} -p {config[macs2_pvalue]} {params.custom} &> {log}
            """

if CONTROL:
    rule call_peaks_Genrich:
        input: control = "04aln_downsample/{control}-downsample.sorted.bam", case="04aln_downsample/{case}-downsample.sorted.bam"
        output: bed = "09peak_Genrich/{case}_vs_{control}_Genrich_peaks.bed"
        log: "00log/{case}_vs_{control}_call_peaks_Genrich.log"
        params:
            name = "{case}_vs_{control}_Genrich",
            jobname = "{case}_vs{control}_Genrich",
            custom = config.get("Genrich_args")
        message: "call_peaks Genrich {input}: {threads} threads"
        shell:
            """
            config[Genrich_path] -t {input.case} \
            -c {input.control} {params.custom} &> {log}
            """
else:
    rule call_peaks_Genrich:
        input: case = "04aln_downsample/{case}-downsample.sorted.bam"
        output: bed = "09peak_Genrich/{case}_Genrich_peaks.bed"
        log: "00log/{case}_call_peaks_Genrich.log"
        params:
            name = "{case}_Genrich",
            jobname = "{case}_Genrich",
            custom = config.get("Genrich_args")
        message: "call_peaks Genrich {input}: {threads} threads"
        shell:
            """
            config[Genrich_path] -t {input.case} \
            {params.custom} &> {log} &> {log}
            """


