from collections import UserDict
from .cytosineread import CytosineRead

class Column(UserDict): 
    def __init__(self,chrom,pos,strand):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand
        super().__init__()

    def __setitem__(self, key, value):
        assert isinstance(value,CytosineRead)
        super().__setitem__(key, value)

    def merge(self, x):
        assert isinstance(x, Column)
        for name,cread in x.items():
            self[name] = cread
