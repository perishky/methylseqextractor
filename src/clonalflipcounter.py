from .windowslider import WindowSlider, WindowIterator


class ClonalFlipCounter:

    def __init__(
        self,
        slider
    ):
        assert isinstance(slider, WindowSlider)
        self.slider = slider

    def iter(self, chrom, start=0, end=None):
        return ClonalFlipCountIterator(self.slider.iter(chrom,start,end))

class ClonalFlipCountIterator:

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
            possible = 0
            flips = 0
            meth = 0
            unmeth = 0
            for clone in pattern["meth"]:
                prev_state = None
                for cpg_state in clone:
                    if cpg_state >= 0:
                        if cpg_state == 1:
                            meth += 1
                        else:
                            unmeth += 1
                        if not prev_state is None:
                            possible += 1
                            if cpg_state != prev_state: flips += 1
                        prev_state = cpg_state
            if possible > 0:
                return { 
                    "chrom": pattern["chrom"],
                    "start": pattern["positions"][0],
                    "end": pattern["positions"][-1],
                    "nsites": len(pattern["positions"]),
                    "nflips":flips, 
                    "npossible": possible,
                    "flip_pct": flips/float(possible),
                    "nmeth": meth,
                    "nunmeth": unmeth,
                    "meth_pct": meth/float(meth+unmeth)
                }
        raise StopIteration
