from .cytosineread import CytosineRead

class Column:
    """
    Provides the reads covering a given CpG site
    """
    def __init__(self,chrom,pos,strand,merge_strands=True):
        """
        attributes:
        - chrom,pos,strand: genomic position of the CpG site
        - merge_strands: whether or not to merge matching CpG sites 
            on opposite strands (default: True)
        """
        self._chrom = chrom
        self._pos = pos
        self._strand = strand
        self._creads = {}
        self.merge_strands = merge_strands

    def get_chrom(self):
        return self._chrom

    def get_pos(self):
        """
        Get the chromosomal position of the CpG site. 
        If the site is on the negative strand and merge_strands is true, 
        then the position is decremented to match the position of the 
        matching CpG site on the positive strand.
        """
        if self.merge_strands and self._strand == "-":
            return self._pos-1
        else:
            return self._pos

    def get_strand(self):
        return self._strand 

    def get_depth(self):
        """ Get the number of reads covering the CpG site """
        return len(self._creads)

    def add_cread(self,cread):
        """ Add a read covering the CpG site """
        assert isinstance(cread,CytosineRead)
        assert self.get_pos() == cread.get_pos()
        self._creads[cread.read.get_name()] = cread

    def get_creads(self):
        """ Get the list of reads covering the CpG site """
        return self._creads.values()

    def merge_negative(self, column):
        """ 
        Assuming the present CpG site is on the positive strand
        and that merge_strands is true, 
        add the reads covering the matching site on the negative strand.
        """
        assert isinstance(column, Column)
        assert column.get_strand() == "-"
        assert self.merge_strands
        assert column.get_pos() == self.get_pos()
        assert self._strand == "+"
        self._strand = "*"
        for cread in column.get_creads():
            self.add_cread(cread)
