from .cpgslider import CpGWindowSlider, CpGIterator
from .windowslider import WindowSlider, WindowIterator
from .window import WindowMaker
from .extractor import Extractor

class CAMDACalculator:

    def __init__(
        self,
        max_fragment,
        extractor,
        merge_strands=True
    ):
        self.max_fragment = max_fragment
        assert isinstance(extractor, Extractor)
        windowmaker = WindowMaker(2*max_fragment,merge_strands)
        windowslider = WindowSlider(windowmaker, extractor)
        self.slider = CpGWindowSlider(windowslider)

    def iter(self, chrom, start=0, end=None):
        iterator = self.slider.iter(chrom,start,end)
        return CAMDAIterator(iterator)

class CAMDAIterator:

    def __init__(
        self,
        iterator
    ):
        assert isinstance(iterator, CpGIterator)
        self.iterator = iterator

    def __iter__(self): 
        return self

    def __next__(self):
        for pattern in self.iterator:
            concurrence_cytosines = 0
            meth_cytosines = 0
            unmeth_cytosines = 0
            meth_clones = 0
            unmeth_clones = 0
            concurrence_clones = 0
            for clone in pattern.meth:
                clone_meth = 0
                clone_unmeth = 0
                for cpg_state in clone:
                    if cpg_state >= 0:
                        if cpg_state == 1:
                            clone_meth += 1
                        else:
                            clone_unmeth += 1
                meth_cytosines += clone_meth
                unmeth_cytosines += clone_unmeth
                if clone_unmeth > 0:
                    if clone_meth > 0:
                        concurrence_cytosines += clone_unmeth
                        concurrence_clones += 1
                    else:
                        unmeth_clones += 1
                elif clone_meth > 0:
                    meth_clones += 1
            number_cytosines = concurrence_cytosines + meth_cytosines + unmeth_cytosines
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
