from .methylseqdataset import MethylSeqDataset
from .window import Window

class ConcurrenceCalculator:

    def __init__(
        self,
        dataset,
        size
    ):
        assert isinstance(dataset,MethylSeqDataset)
        self.window = Window(dataset,size)

    def calculate(self,chrom,start=0,end=None):
        for view in window.slide(chrom,start,end):
            concurrence = 0
            meth = 0
            unmeth = 0
            concurrence_clones = 0
            clones = 0
            for read in view.get_reads():
                clone_meth = 0
                clone_unmeth = 0
                for cread in read.get_creads(view):
                    if cread.is_methylated:
                        clone_meth += 1
                    else:
                        clone_unmeth += 1
                meth += clone_meth
                unmeth += clone_unmeth
                if clone_unmeth > 0:
                    clones += 1
                    if clone_meth > 0:
                        concurrence += clone_unmeth
                        concurrence_clones += 1
                elif clone_meth > 0:
                    clones += 1
            if meth + unmeth > 0:
                yield { 
                    "chrom":view.chrom,
                    "start":view.start,
                    "end":view.end,
                    "nconcurrence":concurrence, 
                    "concurrence": concurrence/float(meth+unmeth),
                    "nconcurrence_clones":concurrence_clones,
                    "unweighted_concurrence":concurrence_clones/float(clones),
                    "nmeth": meth,
                    "meth_pct": float(meth)/float(concurrence+meth+unmeth)
                }
