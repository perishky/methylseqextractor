#!/usr/bin/env python

import os.path
import pandas as pd
import sys
import time
import multiprocessing 
from methylseqextractor import MethylseqExtractor

bamfn = sys.argv[1]
fastafn = sys.argv[2]
outputfn = sys.argv[3]

chromosomes = ["chr"+str(i) for i in range(1,23)] + ["chrX"]

def extract_chromosome_methylation(chrom):
    dat = []
    for cpg in MethylseqExtractor(bamfn, fastafn, chrom):
        dat += cpg
    return pd.DataFrame(dat)

start = time.process_time()

print("Starting extraction ...")

with multiprocessing.Pool(processes=12) as pool:
    for dat in pool.imap(extract_chromosome_methylation, chromosomes):
        if len(dat) > 0:
            dat.to_csv(outputfn + "_" + dat['chrom'][0] + ".csv")

print("Finished in %s seconds." % (time.process_time()-start))
