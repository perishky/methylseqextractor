from .windowslider import WindowSlider, WindowIterator


class ConcurrenceCalculator:

    def __init__(
        self,
        slider
    ):
        assert isinstance(slider, WindowSlider)
        self.slider = slider

    def iter(self, chrom, start=0, end=None):
        slider = self.slider.iter(chrom,start,end)
        return ConcurrenceIterator(slider)

class ConcurrenceIterator:

    def __init__(
        self,
        iterator
    ):
        assert isinstance(iterator, WindowIterator)
        self.iterator = iterator

    def __iter__(self): 
        return self

    def __next__(self):
        for pattern in self.iterator:
            concurrence = 0
            meth = 0
            unmeth = 0
            concurrence_clones = 0
            clones = 0
            for clone in pattern.meth:
                clone_meth = 0
                clone_unmeth = 0
                for cpg_state in clone:
                    if cpg_state >= 0:
                        if cpg_state == 1:
                            clone_meth += 1
                        else:
                            clone_unmeth += 1
                meth += clone_meth
                unmeth += clone_unmeth
                if clone_unmeth > 0:
                    clones += 1
                    if clone_meth > 0:
                        concurrence += clone_unmeth
                        concurrence_clones += 1
                elif clone_meth > 0:
                    clones += 1
            if meth + unmeth > 0:
                return { 
                    "chrom": pattern.chrom,
                    "start": pattern.positions[0],
                    "end": pattern.positions[-1],
                    "nsites": len(pattern.positions),
                    "nconcurrence":concurrence, 
                    "ncytosines":meth+unmeth,
                    "concurrence": concurrence/float(meth+unmeth),
                    "nconcurrence_clones":concurrence_clones,
                    "nclones":clones,
                    "weighted_concurrence":concurrence_clones/float(clones),
                    "nmeth": meth,
                    "nunmeth": unmeth,
                    "meth_pct": float(meth)/float(concurrence+meth+unmeth)
                }
        raise StopIteration
