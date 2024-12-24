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
    def __init__(self,dataset,min_depth=10):
        assert isinstance(dataset,MethylSeqDataset)
        self.dataset = dataset
        self.min_depth = min_depth
        
    def calculate(self,chrom,start=0,end=None):
        iterator = self.dataset.methylation(chrom,start,end)
        for column in iterator:
            meth = 0
            unmeth = 0
            for cread in column.get_creads():
                if cread.is_methylated:
                    meth+=1
                else:
                    unmeth+=1
            if meth+unmeth >= self.min_depth:
                level = float(meth)/float(meth+unmeth)
                yield { 
                    "chrom": column.get_chrom(), 
                    "pos": column.get_pos(),
                    "strand":column.get_strand(),
                    "nmeth":meth,
                    "nunmeth":unmeth,
                    "meth_pct":level
                }
            
