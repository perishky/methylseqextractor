from .methylseqdataset import MethylSeqDataset
from .window import Window

class ConcurrenceCalculator:
    """ 
    Iterate through window views of a given size across the genome
    and calculate the percentage of CpG sites in the view are concurrence CpG sites 
    (unmethylated CpG sites on reads with at least one methylated CpG site 
    also in the view) 
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
          reporting the percentage of CpG sites that are concurrence CpG sites. 
          A dictionary is returned for each view including the following:
          - chrom,start,end: genomic coordinates of the view
          - depth: number of reads in the view
          - nconcurrence: number of concurrence CpG sites in reads
          - concurrence: percentage of CpG sites that are concurrence sites
          - nconcurrence_reads: number of reads with a concurrence CpG site
          - unweighted: percentage of reads that are concurrence reads
          - nmeth: number of methylated CpG sites in the view
          - nunmeth: number of unmethylated CpG sites in the view
          - meth_pct: percentage of CpG sites methylated in the view
        """
        for view in self.window.slide(chrom,start,end):
            concurrence = 0
            meth = 0
            unmeth = 0
            concurrence_clones = 0
            number_clones = 0
            for read in view.get_reads():
                clone_meth = 0
                clone_unmeth = 0
                for cread in read.get_creads(view.start,view.end):
                    if cread.is_methylated:
                        clone_meth += 1
                    else:
                        clone_unmeth += 1
                if clone_meth + clone_unmeth < 2:
                    continue
                meth += clone_meth
                unmeth += clone_unmeth
                if clone_unmeth > 0:
                    number_clones += 1
                    if clone_meth > 0:
                        concurrence += clone_unmeth
                        concurrence_clones += 1
                elif clone_meth > 0:
                    number_clones += 1
            if meth + unmeth > 0:
                yield { 
                    "chrom":view.chrom,
                    "start":view.start,
                    "end":view.end,
                    "depth":number_clones,
                    "nconcurrence":concurrence, 
                    "concurrence": concurrence/float(meth+unmeth),
                    "nconcurrence_reads":concurrence_clones,
                    "unweighted":concurrence_clones/float(number_clones),
                    "nmeth": meth,
                    "nunmeth": unmeth+concurrence,
                    "meth_pct": float(meth)/float(concurrence+meth+unmeth)
                }
