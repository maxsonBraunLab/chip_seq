# common holds pyhton function to be used in the snakefile
import pandas as pd

def gather_files():
    """
    associate samples with fastq files (single or paired) as a dict
    dict = sample:[list of files]
    """
    sample_file_dict = dict()
    all_samples = list( st.index )
    all_files = glob.glob("data/raw/*.fastq.gz")
    for sample in all_samples:
        sample_files = sorted(glob.glob("data/raw/" + sample + "*"))
        sample_file_dict[sample] = sample_files
    return sample_file_dict

def get_control(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """ 
    igg=st.loc[wildcards.sample]['control']
    if igg == "":
        return ""
    else:
        iggbam=f'data/ban/{igg}.ban.sorted.markd.bam'
        return "-c " + iggbam

def get_peak(wildcards):
    """
    return 'broad' flag for macs2 if applicable
    """
    if st.loc[wildcards.sample]['peak'] == "broad":
        return "--broad"
    else:
        return ""

def defect_mode(wildcards, attempt):
    if attempt == 1:
        return ""
    elif attempt > 1:
        return "-D"

def get_callpeaks(wildcards):
    """
    Returns the callpeaks input files
    """
    bam=f"data/markd/{wildcards.sample}.sorted.markd.bam"
    bai=f"data/markd/{wildcards.sample}.sorted.markd.bam.bai"
    if config["USEIGG"]:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/markd/{igg}.sorted.markd.bam'
        iggbam=f'data/markd/{igg}.sorted.markd.bam.bai'
        isigg=config['IGG'] in wildcards.sample
        if not isigg:
            return [bam,bai,iggbam]
        else:
            return [bam,bai]
    else:
        return [bam,bai]
