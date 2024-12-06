from collections import deque 
import pandas as pd

class WindowMaker:
    def __init__(self,size,merge_strands=True):
        self.size = size
        self.merge_strands = merge_strands

    def new_window(self):
        return Window(self.size,self.merge_strands)

class Window:
    def __init__(self,size,merge_strands=True):
        self.size = size
        self.merge_strands = merge_strands
        self.site_reads = deque()
        self.chrom = None
        self.start = None
        self.end = None

    def is_empty(self):
        return len(self.site_reads) == 0

    def wipe(self):
        self.site_reads = deque()

    def add(self,site_read):
        chrom = site_read.get_chrom()
        pos = site_read.get_pos(self.merge_strands)
        if len(self.site_reads) == 0:
            self.chrom = chrom
            self.start = pos
            self.end = pos
            self.site_reads.append(site_read)
            return True
        if chrom == self.chrom \
           and self.end <= pos and pos <= self.start + self.size:
            self.end = pos
            self.site_reads.append(site_read)
            return True
        return False

    def slide(self,site_read):
        if len(self.site_reads) == 0:
            self.add(site_read)
        chrom = site_read.get_chrom()
        pos = site_read.get_pos(self.merge_strands)
        assert chrom != self.chrom or pos >= self.end
        if chrom != self.chrom:
            self.site_reads = deque()
            self.chrom = chrom
            self.start = pos
        self.site_reads.append(site_read)
        self.end = pos
        while pos >= self.start + self.size:
            self.site_reads.popleft()
            self.start = self.site_reads[0].get_pos(self.merge_strands)
        
    def get_pattern(self):
        """
        Returns
        -------
        DNA methylation pattern in the current window view including:
        - chrom: chromosome
        - positions: chromosomal positions of CpG sites 
        - reads: names of reads overlapping the window 
        - meth: matrix of methylation states of CpG sites
          in the window (columns) for each read (rows)
        """
        positions = []
        for site_read in self.site_reads:
            pos = site_read.get_pos(self.merge_strands)
            if len(positions) == 0 or positions[-1] != pos:
                positions.append(pos)
        meth_by_read = dict()
        for site_read in self.site_reads:
            name = site_read.get_name()
            if not name in meth_by_read:
                meth_by_read[name] = [-1]*len(positions)
        idx = 0
        for site_read in self.site_reads:
            name = site_read.get_name()
            pos = site_read.get_pos(self.merge_strands)
            status = 1 if site_read.is_methylated() else 0
            while positions[idx] < pos: idx += 1
            meth_by_read[name][idx] = status
        return Pattern(
            chrom=self.chrom, 
            positions=positions, 
            reads=[name for name in meth_by_read.keys()],
            meth=[x for x in meth_by_read.values()]
        )

class Pattern:
    def __init__(self,chrom,positions,reads,meth):
        assert isinstance(reads,list)
        assert isinstance(meth,list)
        assert len(reads) == len(meth)
        assert len(positions) == len(meth[0])
        self.chrom=chrom
        self.positions=positions
        self.reads=reads
        self.meth=meth

    def __str__(self):
        meth = self.meth
        meth = [["" if m < 0 else str(m) for m in clone] for clone in meth]
        df = pd.DataFrame(meth)
        df.insert(0,"read",self.reads)
        out = ("chrom = " + self.chrom 
                + " positions = " + str(self.positions) + "\n" 
                +  df.to_string(header=False,index=False))
        return out
    
