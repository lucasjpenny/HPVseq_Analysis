#######################
## loading rule modules
include: "rules/common.smk"         ## chekc & install all extra env
include: "rules/genome_coverage.smk"  
include: "rules/end_motif.smk"    
include: "rules/griffin_GC_correction.smk"  
include: "rules/frag_length.smk"   

configfile: "config/config.yaml"

####################
# Read samples files
samples = get_all_samplesnames_list()
print(samples)
###################
## targeted outputs
if config['mappability_correction']: 
    rule all:
        input:
            # expand(config['out_dir'] + "genome_coverage/"+ "fragment_length_no_filter/" + "{sample}.noFilter.fragLen", sample=samples, out_dir=config['out_dir']),
            # expand(config['out_dir']+ "genome_coverage/"+ "depth/" +"{sample}.depth", sample=samples, out_dir=config['out_dir']),
            # expand(config['out_dir']+ "genome_coverage/"+ "collate_and_validate/" + "validate_collate_{sample}", sample=samples, out_dir=config['out_dir']),
            # config['out_dir'] + "end_motif/MDS_scores.txt",
            config['out_dir'] + "fragment_metrics/"+ "fragment_length/" +"fragment_length_summary.csv"
            # Griffin (un)comment out if I want to run or not
            # expand("{out_dir}griffin/GC_mappability/GC_counts/{samples}.GC_counts.txt", samples=samples, out_dir=config['out_dir']), #GC read counts
            # expand("{out_dir}griffin/GC_mappability/GC_bias/{samples}.GC_bias.txt", samples=samples, out_dir=config['out_dir']),
            # expand("{out_dir}griffin/GC_mappability/GC_plots/{samples}.GC_bias.summary.pdf", samples=samples, out_dir=config['out_dir'])
            

#GRIFFIN to delete
# expand("{out_dir}griffin/GC_mappability/mappability_bias/{samples}.mappability_bias.txt", samples=samples, out_dir=config['out_dir']),
# expand("{out_dir}griffin/GC_mappability/mappability_plots/{samples}.mappability_bias.pdf", samples=samples, out_dir=config['out_dir']),
# expand("{out_dir}griffin/GC_mappability/mappability_plots/{samples}.mappability_bias.read_coverage_distribution.pdf", samples=samples, out_dir=config['out_dir']),
# expand("{out_dir}griffin/GC_mappability/samples.GC.yaml", out_dir=config['out_dir'])