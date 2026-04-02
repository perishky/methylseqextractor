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
        positions = deque()
        columns = self.dataset.methylation(chrom,start,end)
        for column in columns:
            positions.append(column.get_pos())
            if column.get_pos() > win_end-1:
                if len(reads) >= self.min_depth:
                    yield WindowView(chrom,win_start,win_end,reads,positions)
                win_end = column.get_pos()+1 
                win_start = win_end - self.size + 1
                while len(positions) > 0:
                    if positions[0] < win_start:
                        positions.popleft()
                    else:
                        break
                while len(reads) > 0:
                    if reads[0].get_end() < win_start+1:
                        reads.popleft()
                    else:
                        break
            for cread in column.get_creads():
                if cread is cread.read.get_leftmost():
                    reads.append(cread.read)
        if len(reads) >= self.min_depth:
            yield WindowView(chrom,win_start,win_end,reads,positions)
        
class WindowView:
    """
    Snapshot view of the DNA methylation statuses of CpGs
    in reads within a given region of the genome
    """
    def __init__(self,chrom,start,end,reads,positions):
        """
        arguments:
        - chrom,start,end: genomic region of interest
        - reads: list of reads in the region (type Read) 
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.reads = reads
        self.positions = positions

    def get_reads(self):
        """
        returns: The list of reads overlapping the genomic region
        """
        return [read for read in self.reads
                if read.get_end()>self.start
                and read.get_start()<self.end]

    def get_grid(self):
        """
        returns: 2D list of DNA methylation statuses of CpGs on the reads 
        overlapping the genomic region. Each row corresponds to a read and each 
        column corresponds to a CpG site. The value is "1" if the CpG is methylated
        in the read, "0" if it is unmethylated or "" if the read does not cover the CpG site.
        """
        grid = []
        for read in self.reads:
            read_meth = [""]*len(self.positions)
            for cread in read.get_creads(self.start,self.end):
                try:
                    idx = self.positions.index(cread.get_pos())
                except ValueError:
                    print("Error: cread pos not in positions")
                    print(self.positions)
                    print(cread.get_pos())
                    continue
                read_meth[idx] = "1" if cread.is_methylated else "0"
            grid += [read_meth]
        return grid
    
    def __str__(self):
        """
        returns: text description of the DNA methylation statuses
        of CpGs on the reads overlapping the genomic region
        """
        meth = pd.DataFrame(self.get_grid())
        meth.insert(0,"read",[read.get_name() for read in self.reads])
        return (
            "chrom = " + self.chrom 
            + ":" + str(self.start) + "-" + str(self.end) 
            + "\npositions = \n " 
            + "\n ".join([str(pos) for pos in self.positions])
            + "\nmeth=\n" 
            +  meth.to_string(header=False,index=False))
        
