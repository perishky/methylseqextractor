from .extractor import Extractor, ExtractionIterator
from .window import WindowMaker, Window
from .siteread import SiteRead

class WindowSlider: 

    """Slides a window across the genome showing methylation patterns"""

    def __init__(self, windowmaker, extractor):
        """
        Parameters
        ----------
        windowmaker: WindowMaker
        extractor: Extractor
        """
        assert isinstance(windowmaker, WindowMaker)
        assert isinstance(extractor, Extractor)
        self.windowmaker = windowmaker
        self.extractor = extractor

    def iter(self, chrom, start=0, end=None):
        """
        Parameters
        ----------
        chrom:
        start:
        end:

        Returns
        -------
        WindowIterator
        """
        return WindowIterator(
            self.windowmaker.new_window(),
            self.extractor.iter(chrom,start,end))

class WindowIterator:
    def __init__(self, window, iterator):
        assert isinstance(window, Window)
        assert isinstance(iterator, ExtractionIterator)
        self.window = window
        self.iterator = iterator

    def __iter__(self):
        return self

    def __next__(self):
        try:
            for site_reads in self.iterator: 
                for site_read in site_reads:
                    if not self.window.add(site_read):
                        pattern = self.window.get_pattern()
                        self.window.slide(site_read)
                        for site_read in site_reads:
                            self.window.add(site_read)
                        return pattern
        except StopIteration:
            pass
        if not self.window.is_empty():
            pattern = self.window.get_pattern()
            self.window.wipe()
            return pattern
        else:
            raise StopIteration
