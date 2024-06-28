from collections import deque 

class MethylationWindow:
    def __init__(self,size,merge_strands=True):
        self.size = size
        self.merge_strands = merge_strands
        self.site_reads = deque()
        self.chrom = None
        self.start = None
        self.stop = None

    def is_empty(self):
        return len(self.site_reads) == 0

    def add(self,site_read):
        chrom = self.site_read.get_chrom()
        pos = self.site_read.get_pos(self.merge_strands)
        if len(self.site_reads) == 0:
            self.chrom = chrom
            self.start = pos
            self.stop = pos
            self.site_reads.appendleft(site_read)
            return True
        if chrom == self.chrom \
           and self.end <= pos and pos <= self.start + self.size:
            self.stop = pos
            self.site_reads.appendleft(site_read)
            return True
        return False

    def shift(self,site_read):
        if len(self.site_reads) == 0:
            self.add(site_read)
        chrom = self.site_read.get_chrom()
        pos = self.site_read.get_pos(self.merge_strands)
        assert chrom != self.chrom or pos >= self.stop
        if chrom != self.chrom:
            self.site_reads = deque()
            self.chrom = chrom
            self.start = pos
        self.site_reads.appendleft(site_read)
        self.stop = pos
        while pos >= self.start + self.size:
            self.site_reads.pop()
            self.start = self.site_reads[0].get_pos(self.merge_strands)
        
    def get_pattern(self):
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
        if site_read in site_reads:
            name = site_read.get_name()
            pos = site_read.get_pos(self.merge_strands)
            status = 1 if site_read.is_methylated() else 0
            while positions[idx] < pos: idx += 1
                idx += 1
            meth_by_read[name][idx] = status
        return {
            "chrom": self.chrom,
            "pos": positions,
            "reads": meth_by_read.keys(),
            "meth": [x for x in meth_by_read.values()]
        }
            

    
