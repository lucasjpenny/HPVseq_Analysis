analysis='hematopoetic'
CPU=1
mem='8G'

griffin= config['griff_dir']
basedir= config['out_dir'] + "griffin/"
ref= config['genome_ref']
outdir=basedir + '/output/nucleosome_profiling/' + analysis

sites=griffin + '/site_configs/' + analysis + '_sites.yaml'
counts=basedir + '/output/GC_correction'
encode_exclude=griffin + '/Ref/encode_unified_GRCh38_exclusion_list.bed'
centromere_path=griffin + '/Ref/hg38_centromeres.bed'
gap_path=griffin + '/Ref/hg38_gaps.bed'
patch_path=griffin + '/Ref/hg38_fix_patches.bed'
alternative_haplotype_path=griffin + '/Ref/hg38_alternative_haplotypes.bed'


rule calc_cov:
	input:
		bam = lambda wildcards: config["samples"][wildcards.samples]['bam'],
		GC_bias = lambda wildcards: config["samples"][wildcards.samples]['GC_bias']
		bam =lambda wildcards: get_all_samples_list(),

	output:
		uncorrected_bw = temp(config['tmp_dir']+"/{samples}/tmp_bigWig/{samples}.uncorrected.bw"),
		GC_corrected_bw = temp(config['tmp_dir']+"/{samples}/tmp_bigWig/{samples}.GC_corrected.bw"),
		tmp_pybedtools = temp(directory(config['tmp_dir']+"/{samples}/tmp_pybedtools"))
	params:
		sample_name = "{samples}",
		mappability_bias = 'none',
		mappability_correction = 'False',

		tmp_dir=config['tmp_dir'],

		reference_genome = config['reference_genome'],
		mappability_bw = config['mappability_bw'],
		chrom_sizes_path = config['chrom_sizes_path'],

		sites_yaml = config['sites_yaml'],
		griffin_scripts_dir = config['griffin_scripts_dir'],
		griffin_coverage_script = config['griffin_scripts_dir']+'/griffin_coverage.py',

		chrom_column=config['chrom_column'],
		position_column=config['position_column'],
		strand_column=config['strand_column'],
		chroms = config['chroms'],

		norm_window = config['norm_window'],
		size_range=config['size_range'],
		map_quality=config['map_quality'],

		number_of_sites=config['number_of_sites'],
		sort_by=config['sort_by'],
		ascending=config['ascending'],

		CPU = config['calc_cov']['ncpus']
		
	shell:
		"time {params.griffin_coverage_script} \
		--sample_name {params.sample_name} \
		--bam {input.bam} \
		--GC_bias {input.GC_bias} \
		--mappability_bias {params.mappability_bias} \
		--mappability_correction {params.mappability_correction} \
		--tmp_dir {params.tmp_dir} \
		--reference_genome {params.reference_genome} \
		--mappability_bw {params.mappability_bw} \
		--chrom_sizes_path {params.chrom_sizes_path} \
		--sites_yaml {params.sites_yaml} \
		--griffin_scripts {params.griffin_scripts_dir} \
		--chrom_column {params.chrom_column} \
		--position_column {params.position_column} \
		--strand_column {params.strand_column} \
		--chroms {params.chroms} \
		--norm_window {params.norm_window} \
		--size_range {params.size_range} \
		--map_quality {params.map_quality} \
		--number_of_sites {params.number_of_sites} \
		--sort_by {params.sort_by} \
		--ascending {params.ascending} \
		--CPU {params.CPU} "