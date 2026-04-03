from .methylseqdataset import MethylSeqDataset
from .window import Window

class ClonalFlipCounter:
    """
    Iterate through window views of the genome of a given size 
    and calculate the percentage of times DNA methylation changes 
    across the reads in each view relative to the number
    of possible changes.
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
          - depth: number of reads in the view
          - num: number of UM flips in the view reads
          - nmu: number of MU flips in the view reads
          - nmm: number of MM pairs in the view reads
          - nuu: number of UU pairs in the view reads
          - nflips: number of flips in the view reads (UM + MU)
          - npossible: number of possible flips in the view reads = 2 x min(nmeth,nunmeth)
          - flip_pct: percentage of flips in the view reads = nflips/npossible
          - nmeth: number of methylated CpG sites in the view
          - nunmeth: number of unmethylated CpG sites in the view
          - meth_pct: percentage of CpG sites methylated in the view
        """
        for view in self.window.slide(chrom,start,end):
            meth = 0
            unmeth = 0
            um_pairs = 0
            mu_pairs = 0
            mm_pairs = 0
            uu_pairs = 0
            reads = view.get_reads()
            for read in reads:
                prev_methylated = None
                for cread in read.get_creads(view.start,view.end):
                    if not prev_methylated is None:
                        if prev_methylated:
                            if cread.is_methylated:
                                mm_pairs += 1
                                meth += 1
                            else:
                                mu_pairs += 1
                                unmeth += 1
                        else:
                            if cread.is_methylated:
                                um_pairs += 1
                                meth += 1
                            else:
                                uu_pairs += 1
                                unmeth += 1
                    prev_methylated = cread.is_methylated
            ## claculate min of int variables meth and unmeth
            possible = min(meth,unmeth)*2
            if possible > 0:
                yield { 
                    "chrom": view.chrom,
                    "start": view.start,
                    "end": view.end,
                    "depth": len(reads),
                    "num": um_pairs,
                    "nmu": mu_pairs,
                    "nmm": mm_pairs,
                    "nuu": uu_pairs,
                    "nflips": um_pairs + mu_pairs,
                    "npossible": possible,
                    "flip_pct": float(um_pairs + mu_pairs)/float(possible),
                    "nmeth":meth,
                    "nunmeth":unmeth,
                    "meth_pct": float(meth)/float(meth+unmeth)
                }
