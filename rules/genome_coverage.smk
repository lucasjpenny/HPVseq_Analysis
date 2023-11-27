# ------------------------------ #
#  Beginning of Snakemake Rules  #
# ------------------------------ #

path_to_data = config['out_dir'] 

rule proper_pair_filter:
    input:
        lambda wildcards: get_sample_path("OPC_test", wildcards.sample)
    output:
        temp_bam=temp(path_to_data + "genome_coverage/"+ "properPaired/" + "{sample}.properPaired.bam"),
        sorted_bam=path_to_data + "genome_coverage/"+ "properPaired/" + "{sample}.properPaired.sorted.bam"
    params:
        mapq=config['mapq'],
        tmppath = path_to_data + "genome_coverage/"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        """
        samtools view -hf 2 -q {params.mapq} {input} | samtools view -bhSo {output.temp_bam}
        samtools sort -T {params.tmppath} -o {output.sorted_bam} {output.temp_bam}
        samtools index {output.sorted_bam}
        """

rule calculate_depth:
    input:
        sorted_bam=path_to_data + "genome_coverage/"+ "properPaired/" + "{sample}.properPaired.sorted.bam"
    output:
        depth=path_to_data + "genome_coverage/"+ "depth/" +"{sample}.depth"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        "samtools depth {input.sorted_bam} > {output.depth}"

rule collate_and_validate:
    input:
        sorted_bam=path_to_data + "genome_coverage/"+ "properPaired/" + "{sample}.properPaired.sorted.bam"
    output:
        collate=path_to_data + "genome_coverage/"+ "collate_and_validate/" + "collate_{sample}",
        validate=path_to_data + "genome_coverage/"+ "collate_and_validate/" + "validate_collate_{sample}"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        """
        samtools collate -o {output.collate} {input.sorted_bam}
        samtools view {output.collate} | awk '{{print $1}}' | sort | uniq -c | awk '{{print $1}}' | sort | uniq -c > {output.validate}
        """

rule fragment_length_no_filter:
    input:
        sorted_bam=lambda wildcards: path_to_data + "genome_coverage/" + "properPaired/" + wildcards.sample + ".properPaired.sorted.bam"
    output:
        frag_len=path_to_data + "genome_coverage/" + "fragment_length_no_filter/" + "{sample}.noFilter.fragLen"
    params:
        tmp1=path_to_data + "genome_coverage/" + "properPaired/" + "{sample}.tmp1",
        tmp2=path_to_data + "genome_coverage/" + "properPaired/" + "{sample}.tmp2",
        tmp=path_to_data + "genome_coverage/" + "properPaired/" + "{sample}.tmp"
    resources: 
        cpus=1, 
        mem_mb=15000, 
        time_min=60
    shell:
        """
        samtools view {input.sorted_bam} | awk 'NR%2==1{{s=$1;for(i=2;i<=12;i++){{s=s"\t"$i}} print s}}' > {params.tmp1}
        samtools view {input.sorted_bam} | awk 'NR%2==0{{s=$1;for(i=2;i<=12;i++){{s=s"\t"$i}} print s}}' > {params.tmp2}
        paste -d '\t' {params.tmp1} {params.tmp2} > {params.tmp}
        echo -e 'start\tend\tfragment_length' > {output.frag_len}
        awk -F'\t' '{{OFS="\t"; if($9>0){{nseq2=length($22); end=$16+nseq2-1; print $4, end, $9}}else{{nseq2=length($10); end=$4+nseq2-1; print $16, end, $21}}}}' {params.tmp} >> {output.frag_len}
        rm {params.tmp1} {params.tmp2} {params.tmp}
        """



