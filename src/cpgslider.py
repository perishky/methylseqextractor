from .windowslider import WindowIterator, WindowSlider
from .window import Pattern

import pandas as pd


class CpGWindowSlider: 

    """
    Slides a window centered at CpG sites across 
    the genome showing methylation patterns in reads 
    that contain the CpG site
    """

    def __init__(self, windowslider):
        assert isinstance(windowslider, WindowSlider)
        self.windowslider = windowslider

    def iter(self, chrom, start=0, end=None):
        """
        Parameters
        ----------
        chrom:
        start:
        end:

        Returns
        -------
        CpGIterator
        """
        return CpGIterator(self.windowslider,chrom,start,end)

class CpGIterator:
    def __init__(self, slider,chrom, start, end):
        assert isinstance(slider, WindowSlider)
        self.size = slider.windowmaker.size
        region_start = start-self.size//2
        region_end = end+self.size//2 if not end is None else None
        self.iterator = slider.iter(chrom,region_start,region_end)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.current = start-1
        self.pattern = None
        self.next_pattern = None
        try:
            self.pattern = next(self.iterator)
            self.next_pattern = next(self.iterator)
        except StopIteration:
            pass
        
    def __iter__(self):
        return self

    def __next__(self):
        if self.pattern is None:
            raise StopIteration
        self.current += 1
        for i in range(len(self.pattern.positions)):
            if self.current < self.pattern.positions[i]:
                self.current = self.pattern.positions[i]
                break
        if not self.end is None and self.current > self.end \
           or self.current > self.pattern.positions[-1]:
            self.pattern = None
            raise StopIteration
        window_start = self.current-self.size//2
        if window_start > self.pattern.positions[0]:
            if not self.next_pattern is None:
                self.pattern = self.next_pattern
                self.next_pattern = None 
                try:
                    self.next_pattern = next(self.iterator)
                except StopIteration:
                    pass
        return CpGPattern(self.pattern,self.current)

class CpGPattern (Pattern):
    def __init__(self,pattern,position):
        assert isinstance(pattern,Pattern)
        assert position in pattern.positions
        self.position = position
        cpg_idx = pattern.positions.index(position)
        ## identify reads that include the CpG site
        reads = []
        meth = []
        for i in range(len(pattern.reads)):
            if pattern.meth[i][cpg_idx] >= 0:
                reads += [pattern.reads[i]]
                meth += [pattern.meth[i]]
        ## identify the first CpG in the retained reads
        found = False
        for first_idx in range(len(pattern.positions)):
            for clone in meth:
                if clone[first_idx] >= 0:
                    found = True
                    break
            if found: 
                break
        ## identify the last CpG in the retained reads
        found = False
        for last_idx in range(len(pattern.positions)-1,-1,-1):
            for clone in meth:
                if clone[last_idx] >= 0:
                    found = True
                    break
            if found:
                break
        assert first_idx <= last_idx
        ## keep info for only the CpG sites in the retained reads
        keep = slice(first_idx,last_idx+1)
        meth = [clone[keep] for clone in meth]
        positions = pattern.positions[keep]
        super().__init__(
            chrom=pattern.chrom,
            positions=positions,
            reads=reads,
            meth=meth
        )
        
    def __str__(self):
        meth = self.meth
        meth = [["" if m < 0 else str(m) for m in clone] for clone in meth]
        meth = pd.DataFrame(meth)
        meth.insert(0,"read",self.reads)
        positions = [str(p) for p in self.positions]
        cpg_idx = positions.index(str(self.position))
        positions[cpg_idx] = "*"+positions[cpg_idx]+"*"
        for i in range(len(meth)):
            meth[i] = [str(m) for m in meth[i]]
            meth[i][cpg_idx] = "*"+meth[i][cpg_idx]+"*"
        out = ("chrom = " + self.chrom 
                + " positions = " + positions + "\n" 
                +  meth.to_string(header=False,index=False))
        return out
