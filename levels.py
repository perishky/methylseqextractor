def calculate(cpg):
    meth = 0
    unmeth = 0
    for read in cpg:
        if read['is_methylated']:
            meth++
        else:
            unmeth++
    return { 
        "chrom":read['chrom'], 
        "pos":read['pos'],
        "meth":meth,
        "unmeth":unmeth,
        "level":float(meth)/float(meth+unmeth)
    }
