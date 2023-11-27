#  Beginning of Snakemake Rules  #
# ------------------------------ #

path_to_data = config['out_dir'] 


import pandas as pd

# Read the samples file
samples = pd.read_csv("samples.csv")

# A function to get the paths of raw files for a given sample
def get_raw_files(sample):
    return [f"/cluster/projects/scottgroup/data/230504_A00827_0793_BHTY2LDSX5_Bratman_Tracy/{sample}_R1_001.fastq.gz",
            f"/cluster/projects/scottgroup/data/230504_A00827_0793_BHTY2LDSX5_Bratman_Tracy/{sample}_R2_001.fastq.gz"]

# Rule for UMI tools extract
rule umitools_extract:
    input:
        lambda wildcards: get_raw_files(wildcards.sample)
    output:
        r1_extracted = "/cluster/projects/scottgroup/people/lucas/HPV-seq/OPC/bwa/../fastq/extracted.{sample}_R1_001.fastq.gz",
        r2_extracted = "/cluster/projects/scottgroup/people/lucas/HPV-seq/OPC/bwa/../fastq/extracted.{sample}_R2_001.fastq.gz"
    resources:
        cpus=1,
        mem_mb=40000,
        time_min=200
    shell:
        "cd /cluster/projects/scottgroup/people/lucas/HPV-seq/OPC/bwa/../fastq/ && "
        "umi_tools extract --extract-method=regex -L {output.r1_extracted}.log --bc-pattern='(?P<umi_1>^[ACGT]{6})' "
        "--filtered-out=/cluster/projects/scottgroup/people/lucas/HPV-seq/OPC/bwa/../fastq/nonUMI.{wildcards.sample}.filtered1.fastq "
        "--filtered-out2=/cluster/projects/scottgroup/people/lucas/HPV-seq/OPC/bwa/../fastq/nonUMI.{wildcards.sample}.filtered2.fastq "
        "--stdin={input[0]} "
        "--stdout={output.r1_extracted} "
        "--read2-in={input[1]} "
        "--read2-out={output.r2_extracted}"





rule process_reads:
    input:
        bam_sorted=lambda wildcards: get_sample_path("OPC_test", wildcards.sample)
    output:
        dedup_log=path_to_data + "fragment_metrics/"+ "virusDuplex/mapped_only/" + "{sample}_dedup.log",
        dedup_bam=temp(path_to_data + "fragment_metrics/"+ "virusDuplex/mapped_only/" + "{sample}_HPVUnmapped_Hg19.Aligned_sorted_dedup.bam"),
        filtered_bam=path_to_data + "fragment_metrics/"+ "virusDuplex/mapped_only/" + "{sample}_HPVUnmapped_Hg19.Aligned_sorted_dedup.filtermapq.bam"
    resources: cpus=1, mem_mb=40000, time_min=200
    shell:
        """
        #module load samtools/1.9
        umi_tools dedup --stdin={input.bam_sorted} --log={output.dedup_log} --unmapped-reads=discard --paired > {output.dedup_bam}
        samtools view -bq 30 {output.dedup_bam} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """

rule align_reads:
    input:
        unmapped_bam=path_to_data + "fragment_metrics/"+ "virusDuplex/mapped_only/" + "{sample}_HPVUnmapped_Hg19.Aligned_sorted_dedup.filtermapq.bam"
    output:
        aligned_bam=path_to_data + "fragment_metrics/virusDuplex/fusion_aligned/{sample}_fusion_aligned.bam"
    params:
        genome_ref=config['genome_ref']
    resources: 
        cpus=8, 
        mem_mb=40000, 
        time_min=120
    shell:
        """
        # Filter mapped reads, convert to FASTQ, align with BWA, and convert back to BAM
        samtools view -b -@ {resources.cpus} -F 4 {input.unmapped_bam} | \
        samtools collate -Ou - | samtools fastq - | \
        bwa mem -M -p -t {resources.cpus} {params.genome_ref} - | \
        samtools view -b -o {output.aligned_bam} -
        """


rule sortindex_reads:
    input:
        # unmapped_bam=path_to_data + "fragment_metrics/"+ "virusDuplex/mapped_only/" + "{sample}_HPVUnmapped_Hg19.Aligned_sorted_dedup.filtermapq.bam"
        bam = path_to_data + "fragment_metrics/"+ "virusDuplex/fusion_aligned/" + "{sample}_fusion_aligned.bam"
    output:
        aligned_bam=path_to_data + "fragment_metrics/"+ "virusDuplex/fusion_aligned/" + "{sample}_fusion_aligned_sorted.bam"
    params:
        genome_ref = config['genome_ref']
    resources: cpus=8, mem_mb=40000, time_min=120
    shell:
        """
        samtools sort -@ {resources.cpus} -o {output.aligned_bam} {input.bam} && samtools index {output.aligned_bam}
        """

