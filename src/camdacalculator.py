from .cpgslider import CpGWindowSlider, CpGIterator
from .windowslider import WindowSlider, WindowIterator
from .window import WindowMaker
from .extractor import Extractor

class CAMDACalculator:

    def __init__(
        self,
        region_size,
        read_length,
        extractor,
        merge_strands=True
    ):
        self.region_size = region_size
        self.read_length = read_length
        assert isinstance(extractor, Extractor)
        windowmaker = WindowMaker(2*read_length,merge_strands)
        windowslider = WindowSlider(windowmaker, extractor)
        self.slider = CpGWindowSlider(windowslider)

    def iter(self, chrom, start=0, end=None):
        iterator = self.slider.iter(chrom,start,end)
        return CAMDAIterator(iterator,self.region_size)

class CAMDAIterator:

    def __init__(
        self,
        iterator,
        region_size
    ):
        assert isinstance(iterator, CpGIterator)
        self.iterator = iterator
        self.region_size = region_size

    def __iter__(self): 
        return self

    def __next__(self):
        for pattern in self.iterator:
            concurrence_cytosines = 0
            number_cytosines = 0
            meth_clones = 0
            unmeth_clones = 0
            concurrence_clones = 0
            region_start = pattern.position - self.region_size//2
            region_end = region_start + self.region_size 
            for clone in pattern.meth:
                clone_meth = 0
                clone_unmeth = 0
                region_unmeth = 0
                for i in range(len(pattern.positions)):
                    cpg_state = clone[i]
                    position = pattern.positions[i]
                    in_region = position >= region_start and position < region_end
                    if cpg_state >= 0:
                        if in_region:
                            number_cytosines += 1
                        if cpg_state == 1:
                            clone_meth += 1
                        else:
                            clone_unmeth += 1
                            if in_region: region_unmeth += 1
                if clone_unmeth > 0:
                    if clone_meth > 0:
                        concurrence_cytosines += region_unmeth
                        concurrence_clones += 1
                    else:
                        unmeth_clones += 1
                elif clone_meth > 0:
                    meth_clones += 1
            number_clones = concurrence_clones + meth_clones + unmeth_clones
            if number_clones > 0:
                return { 
                    "chrom": pattern.chrom,
                    "pos": pattern.position,
                    "ncytosines":number_cytosines,
                    "camda": concurrence_cytosines/float(number_cytosines),
                    "nclones":number_clones,
                    "weighted_camda":concurrence_clones/float(number_clones)
                }
        raise StopIteration
