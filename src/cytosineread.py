from .read import Read

class CytosineRead:
    def __init__(self,read,pos,base,baseq,is_methylated,merge_strands=True):
        assert isinstance(read, Read)
        self.read = read 
        self._pos = pos
        self.base = base
        self.baseq = baseq
        self.is_methylated = is_methylated
        self.merge_strands = merge_strands

    def get_chrom(self):
        return self.read.get_chrom()

    def get_pos(self):
        if self.merge_strands and self.get_strand() == "-":
            return self._pos-1
        else:
            return self._pos

    def get_strand(self):
        return self.read.get_strand()

    def to_dict(self):
        return {
            "read": self.read.get_name(),
            "chrom": self.get_chrom(), 
            "pos": self.get_pos(), 
            "strand": self.get_strand(),
            "base": self.base, 
            "is_methylated": self.is_methylated
        }

    def __str__(self):
        return str(self.to_dict())
