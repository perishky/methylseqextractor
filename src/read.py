class Read:
    def __init__(self,name,chrom,start,end,strand,quality,merge_strands=True):
        self._name = name
        self._chrom = chrom
        self._start = start
        self._end = end
        self._strand = strand
        self._quality = quality
        self._creads = []
        self.merge_strands = merge_strands

    def get_name(self): return self._name
    def get_chrom(self): 
        return self._chrom
    def get_start(self):
        if self.merge_strands and self._strand == "-":
            return self._start-1
        else:
            return self._start
    def get_end(self):
        return self._end
    def get_strand(self):
        return self._strand
    def get_quality(self):
        return self._quality
    def add_cread(self,cread):
        self._creads.append(cread)
    def get_leftmost(self):
        return self._creads[0]
    def get_creads(self,start=None,end=None):
        if start is None and end is None:
            return self._creads
        else:
            creads = self._creads
            if not start is None:
                creads = [cread for cread in creads if cread.get_pos() > start-1]
            if not end is None:
                creads = [cread for cread in creads if cread.get_pos() < end]
            return creads
    def to_dict(self):
        return {
            "name": self.get_name(),
            "chrom": self.get_chrom(), 
            "start": self.get_start(), 
            "end": self.get_end(),
            "strand": self.get_strand()
        }
    def __str__(self):
        return str(self.to_dict())
