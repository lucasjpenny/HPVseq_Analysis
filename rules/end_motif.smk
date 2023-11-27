# ------------------------------ #
#  Beginning of Snakemake Rules  #
# ------------------------------ #

path_to_data = config['out_dir']
genome_ref = config['genome_ref']


# rule mark_duplicates:
#     input:
#         bam = lambda wildcards: get_sample_path("OPC_test", wildcards.sample)
#     output:
#         dedup_bam = temp(path_to_data +"end_motif/{sample}_deduped.bam"),
#         metrics = path_to_data +"end_motif/{sample}_metrics.txt"
#     resources: cpus=1, mem_mb=15000, time_min=60
#     shell:
#         "java -jar $picard_dir/picard.jar MarkDuplicates "
#         "I={input.bam} "
#         "O={output.dedup_bam} "
#         "M={output.metrics} "
#         "TMP_DIR={path_to_data}end_motif "
#         "REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true"

rule sort_bam:
    input:
        dedup_bam = lambda wildcards: get_sample_path("OPC_test", wildcards.sample)
        # dedup_bam = path_to_data + "end_motif/{sample}_deduped.bam"
    output:
        sorted_bam = temp(path_to_data +"end_motif/{sample}_sorted.bam")
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        "samtools sort -n {input.dedup_bam} -o {output.sorted_bam}"

rule samtools_view:
    input:
        sorted_bam = path_to_data +"end_motif/{sample}_sorted.bam"
    output:
        bedpe = temp(path_to_data +"end_motif/{sample}.bedpe")
    params:
        mapq=config['mapq']
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        "samtools view -bf 0x2 -q {params.mapq} {input.sorted_bam} | bedtools bamtobed -i stdin -bedpe > {output.bedpe}"

rule run_rscript_motif_format:
    input:
        bedpe = path_to_data +"end_motif/{sample}.bedpe"
    output:
        fasta_5= temp(path_to_data + "end_motif/" + "{sample}_fasta_5.bed"),
        fasta_3= temp(path_to_data + "end_motif/" + "{sample}_fasta_3.bed")
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        "Rscript scripts/motif_format_bedpe.R "
        "--id {wildcards.sample} "
        "--bedpe {input.bedpe} "
        "--outdir {path_to_data}end_motif"


rule get_fasta:
    input:
        bed3 = path_to_data +"end_motif/{sample}_fasta_3.bed",
        bed5 = path_to_data +"end_motif/{sample}_fasta_5.bed"
    output:
        fasta_out3 = path_to_data +"end_motif/{sample}_fasta_3_annotated.bed",
        fasta_out5 = path_to_data +"end_motif/{sample}_fasta_5_annotated.bed"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        "bedtools getfasta -bedOut -fi {genome_ref} -bed {input.bed3} > {output.fasta_out3};"
        "bedtools getfasta -bedOut -fi {genome_ref} -bed {input.bed5} > {output.fasta_out5}"

rule run_rscript_motif_get_contexts:
    input:
        fasta_5=lambda wildcards: path_to_data + "end_motif/" + wildcards.sample + "_fasta_5_annotated.bed",
        fasta_3=lambda wildcards: path_to_data + "end_motif/" + wildcards.sample + "_fasta_3_annotated.bed"
    output:
        final_output = path_to_data + "end_motif/{sample}_motifs.txt",
        final_output2 = path_to_data + "end_motif/{sample}_raw.txt",
        final_output3 = path_to_data + "end_motif/{sample}_MDS.txt"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        "Rscript scripts/motif_get_contexts.R "
        "--id {wildcards.sample} "
        "--fasta_5 {input.fasta_5} "
        "--fasta_3 {input.fasta_3} "
        "--outdir {path_to_data}end_motif"


rule extract_second_line:
    input:
        expand(path_to_data + "end_motif/{sample}_MDS.txt", sample=get_all_samplesnames_list())
    output:
        path_to_data + "end_motif/MDS_scores.txt" 
    resources: cpus=1, mem_mb=1000, time_min=5
    shell:
        """
        for file in {input}; do
            sample=$(basename $file _MDS.txt)
            value=$(sed -n '2p' $file)
            echo "$sample $value" >> {output}
        done
        """
