from .cytosineread import CytosineRead

class Column: 
    def __init__(self,chrom,pos,strand,merge_strands=True):
        self._chrom = chrom
        self._pos = pos
        self._strand = strand
        self._creads = {}
        self.merge_strands = merge_strands

    def get_chrom(self): 
        return self._chrom

    def get_pos(self):
        if self.merge_strands and self._strand == "-":
            return self._pos-1
        else:
            return self._pos

    def get_strand(self):
        return self._strand 

    def get_depth(self):
        return len(self._creads)

    def add_cread(self,cread):
        assert isinstance(cread,CytosineRead)
        self._creads[cread.read.get_name()] = cread

    def get_creads(self):
        return self._creads.values()

    def merge_negative(self, column):
        assert isinstance(column, Column)
        assert column.get_strand() == "-"
        assert self.merge_strands
        assert column.get_pos() == self.get_pos()
        assert self._strand == "+"
        self._strand = "*"
        for cread in column.get_creads():
            self.add_cread(cread)
