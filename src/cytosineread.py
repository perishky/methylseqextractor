from .read import Read

class CytosineRead:
    def __init__(self,read,chrom,pos,base,baseq,is_methylated):
        assert isinstance(read, Read)
        self.read = read
        self.chrom = chrom
        self.pos = pos
        self.base = base
        self.baseq = baseq
        self.is_methylated = is_methylated
