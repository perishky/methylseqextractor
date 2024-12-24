from .methylseqdataset import MethylSeqDataset
from .window import Window

class ClonalFlipCounter:
    """
    Iterate through window views of the genome of a given size 
    and calculate the percentage of times DNA methylation changes 
    across the reads in each view.
    """
    def __init__(self,dataset,size,min_depth=10):
        """
        attributes: 
        - dataset: type MethylSeqDataset
        - size: size of each window view in base-pairs
        - min_depth: minimum number of reads in any reported window view
        """
        assert isinstance(dataset,MethylSeqDataset)
        self.window = Window(dataset,size,min_depth=min_depth)

    def calculate(self,chrom,start=0,end=None):
        """
        attributes:
        - chrom,start,end: genomic region of interest

        returns: iterator of window views across the genomic region, 
          reporting the percentage of times DNA methylation changes across 
          the reads in the view. A dictionary is returned for each view including 
          the following:
          - chrom,start,end: genomic coordinates of the view
          - nflips: number of times DNA methylation changed across the reads
          - npossible: the maximum number of times DNAm could have changed
          - flip_pct": nflips/npossible
          - nmeth: number of methylated CpG sites in the view
          - nunmeth: number of unmethylated CpG sites in the view
          - meth_pct: percentage of CpG sites methylated in the view
        """
        for view in self.window.slide(chrom,start,end):
            meth = 0
            unmeth = 0
            possible = 0
            flips = 0
            for read in view.get_reads():
                prev_methylated = None
                for cread in read.get_creads(view.start,view.end):
                    if cread.is_methylated: 
                        meth += 1
                    else:
                        unmeth += 1
                    if not prev_methylated is None:
                        possible += 1
                        if prev_methylated != cread.is_methylated:
                            flips += 1
                    prev_methylated = cread.is_methylated
            if possible > 0:
                yield { 
                    "chrom": view.chrom,
                    "start": view.start,
                    "end": view.end,
                    "nflips":flips, 
                    "npossible":possible,
                    "flip_pct": flips/float(possible),
                    "nmeth":meth,
                    "nunmeth":unmeth,
                    "meth_pct": float(meth)/float(meth+unmeth)
                }
