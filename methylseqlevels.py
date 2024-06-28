from .methylseqextractor import MethylSeqExtractor
from .siteread import SiteRead

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
        return self

    def __next__(self):
        site_reads = next(self.extractor)
        meth = 0
        unmeth = 0
        for site_read in site_reads:
            if site_read.is_methylated():
                meth+=1
            else:
                unmeth+=1
        level = None
        if meth+unmeth > 0:
            level = float(meth)/float(meth+unmeth)
        return { 
            "chrom":site_read.get_chrom(), 
            "pos":site_read.get_pos(),
            "meth":meth,
            "unmeth":unmeth,
            "level":level
        }
