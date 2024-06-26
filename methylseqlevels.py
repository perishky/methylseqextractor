from methylseqextractor import MethylseqExtractor

class MethylSeqLevels: 

    """Iterates through CpG sites and calculates DNAm levels"""

    def __init__(self, extractor):
        """
        Parameter
        ---------
        extractor: MethylSeqExtractor

        Returns
        -------
        Dictionary for the current CpG site including:
        - chromosome ('chromosome')
        - chromosomal position ('pos')
        - number of methylated cytosines ('meth')
        - number of unmethylated cytosines ('unmeth')
        - methylation level ('level')
        """
        self.extractor = extractor

    def __iter__(self):
        return self:

    def __next__(self):
        site = self.extractor.next()        
        meth = 0
        unmeth = 0
        for read in site:
            if read['is_methylated']:
                meth++
            else:
                unmeth++
        level = None
        if meth+unmeth > 0:
            level = float(meth)/float(meth+unmeth)
        return { 
            "chrom":read['chrom'], 
            "pos":read['pos'],
            "meth":meth,
            "unmeth":unmeth,
            "level":level
        }
