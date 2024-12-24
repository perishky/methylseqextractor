from collections import deque 
import pandas as pd

from .methylseqdataset import MethylSeqDataset

class Window:
    """
    Sliding window providing view of the DNA methylation statuses
    of CpG sites in reads across a specified genomic region
    """
    def __init__(self,dataset,size,min_depth=10):
        """
        arguments:
        - dataset: type MethylSeqDataset
        - size: size of window view (basepairs)
        - min_depth: minimum number of reads in any reported window view
        """
        assert isinstance(dataset,MethylSeqDataset)
        self.dataset = dataset
        self.size = size
        self.min_depth = min_depth

    def slide(self,chrom,start=0,end=None):
        """
        Slide the window across the 

        arguments:
        - chrom,start,end: genomic region across which to slide the window

        returns: iterator of window views (type WindowView)
        """
        win_start = start
        win_end = start + self.size
        reads = deque()
        columns = self.dataset.methylation(chrom,start,end)
        for column in columns:
            if column.get_pos() > win_end-1:
                if len(reads) >= self.min_depth:
                    yield WindowView(chrom,win_start,win_end,reads)
                win_end = column.get_pos()+1 
                win_start = win_end - self.size + 1
                while len(reads) > 0:
                    if reads[0].get_end() < win_start+1:
                        reads.popleft()
                    else:
                        break
            for cread in column.get_creads():
                if cread is cread.read.get_leftmost():
                    reads.append(cread.read)
        if len(reads) >= self.min_depth:
            yield WindowView(chrom,win_start,win_end,reads)
        
class WindowView:
    """
    Snapshot view of the DNA methylation statuses of CpGs
    in reads within a given region of the genome
    """
    def __init__(self,chrom,start,end,reads):
        """
        arguments:
        - chrom,start,end: genomic region of interest
        - reads: list of reads in the region (type Read) 
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.reads = reads

    def get_reads(self):
        """
        returns: The list of reads overlapping the genomic region
        """
        return [read for read in self.reads
                if read.get_end()>self.start
                and read.get_start()<self.end]

    def __str__(self):
        """
        returns: text description of the DNA methylation statuses
        of CpGs on the reads overlapping the genomic region
        """
        reads = self.get_reads()
        meth = []
        positions = set()
        for read in reads:
            for cread in read.get_creads(self.start,self.end):
                positions.add(cread.get_pos())
        positions = list(positions)
        positions.sort()
        for read in reads:
            read_meth = [""]*len(positions)
            for cread in read.get_creads(self.start,self.end):
                idx = positions.index(cread.get_pos())
                read_meth[idx] = "1" if cread.is_methylated else "0"
            meth += [read_meth]
        meth = pd.DataFrame(meth)
        meth.insert(0,"read",[read.get_name() for read in reads])
        return (
            "chrom = " + self.chrom 
            + ":" + str(self.start) + "-" + str(self.end) 
            + "\npositions = \n " 
            + "\n ".join([str(pos) for pos in positions])
            + "\nmeth=\n" 
            +  meth.to_string(header=False,index=False))
        
