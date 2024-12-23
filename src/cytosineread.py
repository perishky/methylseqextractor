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
        self.read.creads.append(self)

    def get_dict(self):
        return {
            "read": self.read.name,
            "chrom": self.chrom, 
            "pos": self.pos, 
            "strand": self.read.strand,
            "base": self.base, 
            "is_methylated": self.is_methylated
        }

    def __str__(self):
        return str(self.get_dict())
