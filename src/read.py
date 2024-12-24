class Read:
    """
    Represents an aligned sequencing read providing its genomic coordinates,
    and the DNA methylation statuses of all CpG sites on the read.
    """
    def __init__(self,name,chrom,start,end,strand,quality,merge_strands=True):
        """
        attributes:
        - name: unique Illumina identifier of the aligned read from the BAM file
        - chrom,start,end,strand: genomic coordinates of the aligned read
        - quality: read alignment quality from the BAM file
        - merge_strands: whether or not to merge matching CpG sites 
            on opposite strands (default: True)
        """
        self._name = name
        self._chrom = chrom
        self._start = start
        self._end = end
        self._strand = strand
        self._quality = quality
        self._creads = []
        self.merge_strands = merge_strands

    def get_name(self):
        return self._name
    def get_chrom(self): 
        return self._chrom
    def get_start(self):
        """ 
        Get the genomic start position of the read.
        If the read was on the negative strand and 
        merge_strands is true, then the start position 
        is one less similarly to the positions of cytosines 
        in matching CpG sites on the positive strand. 
        """ 
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
        """
        Add the DNA methylation status of a CpG to the read (type CytosineRead). 
        Assumes that statuses are added in order of genomic position.
        """
        assert (len(self._creads) == 0
                or self._creads[-1].get_pos() < cread.get_pos())
        self._creads.append(cread)
    def get_leftmost(self):
        """ Get the DNA methylation status CpG site with the leftmost position """
        return self._creads[0]
    def get_creads(self,start=None,end=None):
        """
        Get the DNA methylation statuses of CpG sites on the read

        attributes:
        - start,end: restrict chromosomal region of the CpG sites

        returns: List of DNA methylation statuses (type CytosineRead)
        """
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
