configfile: 'config.yaml'
localrules: all

SAMPLES = config['samples']
RAWDATA = config['rawdata']
RESULTS = config['results']
REFERENCE = config['reference']

BAM = expand(f'{RESULTS}/bam-RNA/{{sample}}.bam \
{RESULTS}/bam-RNA/{{sample}}.bam.bai \
{RESULTS}/bam-RNA/samtools-stats/{{sample}}.bam.bc \
{RESULTS}/bam-RNA/samtools-stats/{{sample}}.bam.bc_plots/ \
{RESULTS}/bam-RNA/fastqc-output/{{sample}}_fastqc.html'.split(), sample = SAMPLES)

rule all:
    input:
        BAM

rule hisat2:
    input:
        r1 = f'{RAWDATA}/{{sample}}_R1.fastq',
        r2 = f'{RAWDATA}/{{sample}}_R2.fastq'
    output:
        sam = temp(f'{RESULTS}/hisat2-output/{{sample}}.sam'),
        log = f'{RESULTS}/hisat2-output/.hisat2.{{sample}}.log'
    shell:
        f"""
        hisat2 -p 24 -t --no-unal -x {REFERENCE} -1 {{input.r1}} -2 {{input.r2}} >{{output.sam}} 2>{{output.log}}
        """

rule toBAM:
    input:
        f'{RESULTS}/hisat2-output/{{sample}}.sam'
    output:
        temp(f'{RESULTS}/bam-RNA/{{sample}}.unsorted.bam')
    shell:
        """
        grep -v 'NH:i:[2-9]' {input} | samtools view -h -b -F 4 -q 20 >{output}
        """

rule sortBAM:
    input:
        f'{RESULTS}/bam-RNA/{{sample}}.unsorted.bam'
    output:
        f'{RESULTS}/bam-RNA/{{sample}}.bam'
    shell:
        """
        samtools sort -@ 10 {input} >{output}
        """

rule indexBAM:
    input:
        f'{RESULTS}/bam-RNA/{{sample}}.bam'
    output:
        f'{RESULTS}/bam-RNA/{{sample}}.bam.bai'
    shell:
        """
        samtools index {input}
        """

rule statsBAM:
    input:
        f'{RESULTS}/bam-RNA/{{sample}}.bam'
    output:
        f'{RESULTS}/bam-RNA/samtools-stats/{{sample}}.bam.bc'
    shell:
        """
        samtools stats {input} >{output}
        """

rule plotBAM:
    input:
        f'{RESULTS}/bam-RNA/samtools-stats/{{sample}}.bam.bc'
    output:
        directory(f'{RESULTS}/bam-RNA/samtools-stats/{{sample}}.bam.bc_plots')
    shell:
        """
        plot-bamstats -p {output}/ {input}
        """

rule fastQC:
    input:
        f'{RESULTS}/bam-RNA/{{sample}}.bam'
    output:
        f'{RESULTS}/bam-RNA/fastqc-output/{{sample}}_fastqc.html'
    shell:
        f"""
        fastqc -q -t 24 -o {RESULTS}/bam-RNA/fastqc-output {{input}}
        """
