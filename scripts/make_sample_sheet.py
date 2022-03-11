#!/usr/bin/env python3

import pandas as pd
import sys
import os
import re
import numpy as np

def main():

    # get list of fastq files in raw dir
    raw = sys.argv[1]
    try:
        control = sys.argv[2]
    except IndexError:
        control = "IgG"
    outfile = "config/samplesheet.tsv"

    samps = [os.path.join(raw,f) for f in os.listdir(raw)]
    regex = re.compile(".+\.fastq\.gz")

    # define samples
    sfilt = [s for s in samps if regex.search(s)]
    sfilt = [os.path.basename(x).split(".")[0] for x in sfilt] # rm .fastq.gz suffix
    sfilt = [x.replace("_R1", "").replace("_R2", "") for x in sfilt] # rm R1 R2 info
    sfilt = [x for x in sfilt if control not in x] # rm IgG samples from samples
    sfilt = sorted(sfilt)

    # make df
    df = pd.DataFrame(sfilt, columns = ["sample"])
    df = df.drop_duplicates()
    df["control"] = np.repeat("control_sample", df.shape[0])
    df["peak"] = np.repeat("narrow", df.shape[0])

    # write sample sheet to file
    df.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__": 
    main()
