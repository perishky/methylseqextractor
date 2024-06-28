from .methylseqextractor import MethylSeqExtractor
from .methylationwindow import MethylationWindow
from .siteread import SiteRead

class MethylSeqSlidingWindow: 

    """Slides a window across the genome showing methylation patterns"""

    def __init__(self, window, extractor):
        """
        Parameter
        ---------
        window: MethylationWindow
        extractor: MethylSeqExtractor

        Returns
        -------
        Dictionary for the current window including:
        - chrom: chromosome
        - pos: chromosomal positions of CpG sites 
        - reads: names of reads overlapping the window 
        - meth: matrix of methylation states of CpG sites
          in the window (columns) for each read (rows)
        """
        self.window = window
        self.extractor = extractor
        self.merge_strands = merge_strands

    def __iter__(self):
        return self

    def __next__(self):
        site_reads = []
        for site_reads in self.extractor: 
            added = False
            for site_read in site_reads:
                added = self.window.add(site_read)
                if not added: break
            if not added: break
        if self.window.is_empty(): raise StopIteration
        pattern = self.window.get_pattern()
        for site_read in site_reads:
            self.window.shift(site_read)
        return pattern

