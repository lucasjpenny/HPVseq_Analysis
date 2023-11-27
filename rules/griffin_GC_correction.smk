
griffin= config['griff_dir']
basedir= config['out_dir'] + "griffin/"
ref= config['genome_ref']
outdir=basedir + "GC_mappability/"


rule griffin_mappability_correction:
    input:
        bam_files = lambda wildcards: get_all_samples_list()
    output:
        mappability_bias = outdir + "mappability_bias/" +"{sample}.mappability_bias.txt",
        mappability_plot = outdir + "mappability_plots/" +"{sample}.mappability_bias.read_coverage_distribution.pdf",
        tmp_dir = temp(directory(outdir + "tmp_{sample}/"))
    params:
        griffin = config['griff_dir'],
        outdir = config['out_dir'] + "GC_correction",
        mappability = config['griff_dir'] + "/Ref/k50.Umap.MultiTrackMappability.hg38.bw",
        exclude_paths = config['griff_dir'] + "/Ref/encode_unified_GRCh38_exclusion_list.bed",
        chrom_sizes = config['griff_dir'] + "/Ref/hg38.standard.chrom.sizes",
        map_quality = 20,
        CPU = 8
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        """
        {params.griffin}/scripts/griffin_mappability_correction.py \
        --bam_file {input.bam_files} \
        --bam_file_name {wildcards.sample} \
        --output {output.mappability_bias} \
        --output_plot {output.mappability_plot} \
        --mappability {params.mappability} \
        --exclude_paths {params.exclude_paths} \
        --chrom_sizes {params.chrom_sizes} \
        --map_quality {params.map_quality} \
        --CPU {params.CPU} \
        --tmp_dir {output.tmp_dir}
        """


rule griffin_gc_counts:
    input:
        bam_file = lambda wildcards: get_all_samples_list()
    output:
        gc_counts = outdir + "GC_counts/{sample}.GC_counts.txt"
    params:
        griffin = config['griff_dir'],
        mappable_regions_path = config['griff_dir'] + "/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed",
        ref_seq = config['genome_ref'],
        chrom_sizes = config['griff_dir'] + "/Ref/hg38.standard.chrom.sizes",
        out_dir = outdir,
        map_q = 20,
        size_range = "15 500",
        CPU = 8
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        """
        {params.griffin}/scripts/griffin_GC_counts.py \
        --bam_file {input.bam_file} \
        --bam_file_name {wildcards.sample} \
        --mappable_regions_path {params.mappable_regions_path} \
        --ref_seq {params.ref_seq} \
        --chrom_sizes {params.chrom_sizes} \
        --out_dir {params.out_dir} \
        --map_q {params.map_q} \
        --size_range {params.size_range} \
        --CPU {params.CPU}
        """



rule griffin_gc_bias:
    input:
        bam = outdir + "GC_counts/{sample}.GC_counts.txt"
    output:
        gc_bias = outdir + "GC_bias/{sample}.GC_bias.txt",
        GC_plots_file = outdir + "GC_plots/{sample}.GC_bias.summary.pdf"
    params:
        griffin = config['griff_dir'],
        mappable_name = "k100_minus_exclusion_lists.mappable_regions.hg38",
        genome_GC_frequency = config['griff_dir'] + "/Ref/genome_GC_frequency",
        out_dir = outdir,
        size_range = "15 500"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        """
        {params.griffin}/scripts/griffin_GC_bias.py \
        --bam_file_name {wildcards.sample} \
        --mappable_name {params.mappable_name} \
        --genome_GC_frequency {params.genome_GC_frequency} \
        --out_dir {params.out_dir} \
        --size_range {params.size_range}
        """




