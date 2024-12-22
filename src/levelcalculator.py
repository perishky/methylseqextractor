from .methylseqdataset import MethylSeqDataset
from .cytosineread import CytosineRead
from .read import Read

class LevelCalculator:
    """
    Iterates through CpG sites and calculates DNAm levels

    Parameter
    ---------
    iterator: CytosineIterator
    
    Returns
    -------
    Dictionary for the current CpG site including:
    - chromosome ('chromosome')
    - chromosomal position ('pos')
    - number of methylated cytosines ('meth')
    - number of unmethylated cytosines ('unmeth')
    - methylation level ('level')
    """
    def __init__(self,dataset):
        assert isinstance(dataset,MethylSeqDataset)
        self.dataset = dataset
        
    def calculate(self,chrom,start=0,end=None):
        iterator = self.dataset.methylation(chrom,start,end)
        for column in iterator:
            meth = 0
            unmeth = 0
            for cread in column.values():
                if cread.is_methylated:
                    meth+=1
                else:
                    unmeth+=1
            if meth+unmeth > 0:
                level = float(meth)/float(meth+unmeth)
                strand = cread.read.strand
                if self.dataset.merge_strands: strand = "*"
                yield { 
                    "chrom":cread.read.chrom, 
                    "pos": cread.pos,
                    "strand":strand,
                    "nmeth":meth,
                    "nunmeth":unmeth,
                    "meth_pct":level
                }
            
