from .methylseqdataset import MethylSeqDataset
from .cytosineread import CytosineRead
from .read import Read

class LevelCalculator:
    """
    Iterate through CpG sites and calculate their DNAm levels
    """
    def __init__(self,dataset,min_depth=10):
        """
        attributes:
        - dataset: type MethylSeqDataset
        - min_depth: minimum number of reads covering each reported CpG site
        """
        assert isinstance(dataset,MethylSeqDataset)
        self.dataset = dataset
        self.min_depth = min_depth
        
    def calculate(self,chrom,start=0,end=None):
        """
        attributes:
        - chrom,start,end: genomic region of interest

        returns: iterator of CpG sites across the region 
          and their DNA methylation levels. Sites covered by 
          fewer than min_depth (see init function) are omitted.
          A dictionary is returned for each CpG site including the following:
          - chrom,pos,strand: genomic position of the CpG site
          - nmeth: number of reads in which it is methylated
          - nunmeth: number of reads in which it is unmethylated
          - meth_pct: percentage of reads in which it is methylated
        """
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
            
