from collections import deque 
import pandas as pd

from .methylseqdataset import MethylSeqDataset

class Window:
    def __init__(self,dataset,size):
        assert isinstance(dataset,MethylSeqDataset)
        self.dataset = dataset
        self.size = size

    def slide(self,chrom,start=0,end=None):
        reads = deque()
        iterator = self.dataset.methylation(chrom,start,end)
        win_start = start
        win_end = start + self.size
        for column in iterator:
            if column.pos > win_end:
                if len(reads) > 0:
                    yield WindowView(chrom,win_start,win_end,reads)
                win_end = column.pos 
                win_start = win_end - self.size
                while len(reads) > 0:
                    if reads[0].end < win_start:
                        reads.popleft()
                    else:
                        break
            for cread in column.values():
                if cread.read.creads[0].pos >= column.pos: 
                    reads.append(cread.read)
        if len(reads) > 0:
            yield WindowView(chrom,win_start,win_end,reads)
        
class WindowView:

    def __init__(self,chrom,start,end,reads):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.reads = reads

    def get_reads(self): 
        return [read for read in self.reads \
                if read.end >= self.start and read.start <= self.end]

    def __str__(self):
        reads = self.get_reads()
        meth = []
        positions = set()
        for read in reads:
            for cread in read.get_creads(self.start,self.end):
                positions.add(cread.pos)
        positions = list(positions)
        positions.sort()
        for read in reads:
            read_meth = [""]*len(positions)
            for cread in read.get_creads(self.start,self.end):
                idx = positions.index(cread.pos)
                read_meth[idx] = "1" if cread.is_methylated else "0"
            meth += [read_meth]
        meth = pd.DataFrame(meth)
        meth.insert(0,"read",[read.name for read in reads])
        return (
            "positions = \n " 
            + "\n ".join([str(pos) for pos in positions])
            + "\nmeth=\n" 
            +  meth.to_string(header=False,index=False))
        
