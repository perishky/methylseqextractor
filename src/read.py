from .window import WindowView
from .cytosineread import CytosineRead

class Read:
    def __init__(self,name,chrom,start,end,strand,quality):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.quality = quality
        self.creads = []
    def get_creads(self,view=None):
        if view is None:
            return self.creads
        else:
            assert isinstance(view,WindowView)
            return [cread for cread in self.creads 
                    if cread.pos >= view.start and cread.pos <= view.end]

