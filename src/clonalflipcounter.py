from .methylseqdataset import MethylSeqDataset
from .window import Window

class ClonalFlipCounter:

    def __init__(self,dataset,size,min_depth=10):
        assert isinstance(dataset,MethylSeqDataset)
        self.window = Window(dataset,size,min_depth=min_depth)

    def calculate(self,chrom,start=0,end=None):
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
