class Read:
    def __init__(self,name,chrom,start,end,strand,quality):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.quality = quality
        self.creads = []
    def get_creads(self,start=None,end=None):
        if start is None and end is None:
            return self.creads
        else:
            creads = self.creads
            if not start is None:
                creads = [cread for cread in creads if cread.pos > start-1]
            if not end is None:
                creads = [cread for cread in creads if cread.pos < end]
            return creads
    def get_dict(self):
        return {
            "name": self.name,
            "chrom": self.chrom, 
            "start": self.start, 
            "end": self.end,
            "strand": self.strand
        }
    def __str__(self):
        return str(self.get_dict())
