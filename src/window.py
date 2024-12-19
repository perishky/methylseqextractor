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
        positions = deque()
        iterator = self.dataset.methylation(chrom,start,end)
        snap_start = start
        snap_end = start + self.size
        for column in iterator:
            pos = column[0].pos
            if pos > snap_end:
                if len(positions) > 0:
                    yield WindowView(chrom,snap_start,snap_end,positions,reads)
                snap_end = pos
                snap_start = snap_end - self.size
                while len(positions) > 0:
                    if positions[0] < snap_start:
                        positions.popleft()
                    else:
                        break
                while len(reads) > 0:
                    if reads[0].end < snap_start:
                        reads.popleft()
                    else:
                        break
            positions.append(pos)
            for cread in column:
                if cread.read.creads[0].pos == pos: 
                    reads.append(cread.read)
        if len(positions) > 0:
            yield WindowView(chrom,snap_start,snap_end,positions,reads)
        
class WindowView:

    def __init__(self,chrom,start,end,positions,reads):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.positions = positions
        self.reads = reads

    def get_reads(self): 
        return [read for read in self.reads \
                if read.end >= self.start and read.start <= self.end]

    def get_positions(self):
        return [pos for pos in self.positions \
                if pos >= self.start and pos <= self.end]

    def __str__(self):
        positions = self.get_positions()
        reads = self.get_reads()
        meth = []
        for read in reads:
            read_meth = [""]*len(positions)
            for cread in read.get_creads(self.start,self.end):
                read_meth[positions.index(cread.pos)] = "1" if cread.is_methylated else "0"
            meth += [read_meth]
        meth = pd.DataFrame(meth)
        meth.insert(0,"read",[read.name for read in reads])
        return ("positions = \n " + "\n ".join([str(pos) for pos in positions])
                + "\nmeth=\n" +  meth.to_string(header=False,index=False))
        
