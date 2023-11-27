
# ------------ #
#  Base Setup  #
# ------------ #
import pandas as pd

configfile: "config/config.yaml"

path_to_data = config['data']['base_path']


#############################################
## get taget outputs based on the config file
## either for individual samples or aggregate
## all samples listed in sample_aggr.tsv !!!
#############################################

#FIX THIS TO ACCOMODATE WHAT I NEED

def get_rule_all_input():
    ## ensure extra env installed
    extra_env = "extra_env/all_extra_env_installed",

    ## fixed outputs
    #meth_qc = "aggregated/meth_qc.txt",
    aggr_qc = "aggregated/aggr_qc_report.html",
    meta_quant = "aggregated/meth_count.txt.gz",
    meth_filt = "autos_bfilt/meth_count_autos_bfilt.txt.gz",

    ######################################
    ## aggregated outputs for SAMPLES_aggr
    ## paired-end and spike-in
    if config["aggregate"] and config["paired-end"] and config["spike_in"] and config["frag_profile"]:
        #mult_qc = "aggregated/QC_pe/multiqc_report.html",

        ## spike-ins
        spikein_mult_qc = "aggregated_spikein/QC_pe/multiqc_report.html",
        spikein_meth_qc = "aggregated_spikein/meth_qc.txt",
        spikein_meta_quant = "aggregated_spikein/meth_count.txt.gz",

        #fragment profiles
        fp_gc = "aggregated/fragment_profile_GC_corrected_1mb.tsv",        ## GC corrected fragment profile

        return  extra_env + aggr_qc + meta_quant + meth_filt + spikein_mult_qc + spikein_meth_qc + spikein_meta_quant + fp_gc






def get_cohort_data(cohort):
    """Parses the samplesheet for a specific cohort.
    """
    
    samplesheet = pd.read_csv(config['data']['cohorts'][cohort]['samplesheet'], comment='#').drop_duplicates()
    return samplesheet

def get_all_samples(cohort=None):
    """Retrieves all samples to be processed.
    Does so by calling get_cohort_data, and therefore filters out excluded_cases.
    Keyword arguments:
        cohort -- Name of a cohort, OPTIONAL. If not specified, returns all samples
                  across all cohorts.
    """
    
    all_samples = pd.concat([
        get_cohort_data(cohort_name).assign(cohort_name = cohort_name)
        for cohort_name
        in config['data']['cohorts']
        if config['data']['cohorts'][cohort_name]['active']
        ])

    if cohort is None:
        return(all_samples)
    else:
        return(all_samples[all_samples.cohort_name == cohort])

def get_all_samples_list(cohort=None):
    """Returns all samples in list format.
    By calling get_all_samples()
    """
    
    return get_all_samples(cohort).path.unique().tolist() #CHANGED TO HAVE PATH

def get_all_samplesnames_list(cohort=None):
    """Returns all samples in list format.
    By calling get_all_samples()
    """
    
    return get_all_samples(cohort).sample_name.unique().tolist() #CHANGED To NAMES



def get_all_samples_with_path():
    """Returns a list of tuples with cohort name and sample name."""
    samples = get_all_samples()#[['sample_name','path']]
    return(zip(
        samples.drop_duplicates().sample_name.tolist(),
        samples.drop_duplicates().path.tolist()
    ))

def get_sample_path(cohort, sample):
    """Retrieves the path to the fastq file.
    Keyword arguments:
        cohort -- name of the cohort whose samplesheet should be accessed.
        sample -- identifier of the sample as specified in the samplesheet.
        library -- integer representing the library index as specified in the samplesheet.
        read_in_pair -- 1 or 2 - representing read 1 or read 2 in paired end data.
    """

    cohort_data = get_cohort_data(cohort)
    sample_line = cohort_data[
        (cohort_data.sample_name == sample)]
    return sample_line.path.to_list()[0]

def extract_bedfile_names(directory):
    bed_files = []
    for filename in os.listdir(directory):
        if filename.endswith(".bed"):
            bed_files.append(filename)
    return bed_files