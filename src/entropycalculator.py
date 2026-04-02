from .methylseqdataset import MethylSeqDataset
from .window import Window
import numpy as np

class EntropyCalculator:
    """
    Iterate through window views of the genome of a given size 
    and calculate the DNA methylation similarity-sensitive entropy
    in each view.

    WARNING: This is very slow, e.g. runtime for one sample was 2 hours
    for 100bp windows. Compare this the ClonalFlipCounter which took
    12 minutes for the same sample and window size.

    WARNING: I'm also not sure how informative this measure is beyond
    what could be obtained from methylation levels, i.e. entropy is
    naturally higher when methylation levels are around 50% and lowest
    around 0% and 100%. If we removed the log from the calculation, then
    the resulting measure is just the mean of the similarities between
    all pairs of measurements for each CpG site. This can be calculated
    by just knowing the numbers of methylated and unmethylated measurements
    for each CpG site, i.e. read level information is not needed. Unclear
    how best to increase the contribution of reads to the calculation.
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
          reporting the similarity-sensitive entropy across the reads in the view.
          A dictionary is returned for each view including the following:
          - chrom,start,end: genomic coordinates of the view
          - entropy: similarity-sensitive entropy of the DNAm in the view
          - depth: number of reads in the view
          - sites: number of CpG sites in the view
          - comps: number of pairwise site comparisons made between reads in the view
        """
        for view in self.window.slide(chrom,start,end):
            grid = view.get_grid()
            mean_sim = [0]*len(grid)
            ncomps = 0
            for i in range(len(grid)):
                sim = []
                for j in range(len(grid)):
                    n = 0
                    d = 0
                    for k in range(len(grid[i])):
                        if grid[i][k] != "" and grid[j][k] != "":
                            n += 1
                            ncomps += 1
                            if grid[i][k] != grid[j][k]:
                                d += 1
                    if n > 0:
                        sim.append(1-(float(d)/float(n)))
                if len(sim) > 0:
                    mean_sim[i] = np.mean(sim)
            entropy = -np.mean([np.log2(np.clip(x, 1e-10,1.0)) for x in mean_sim])
            yield { 
                "chrom": view.chrom,
                "start": view.start,
                "end": view.end,
                "entropy": entropy,
                "depth": len(grid),
                "sites": len(grid[0]),
                "comps": ncomps
            }


