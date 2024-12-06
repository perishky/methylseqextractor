import pysam
import numpy
from .siteread import SiteRead

class Extractor:

    """Iterates through methyl-seq reads by CpG site"""

    def __init__(
        self, 
        bamfn,
        fastafn,
        min_mapq=10, 
        min_read_quality=5,
        min_base_quality=5,
        min_depth=10
    ):
        """ 
        Parameters
        ----------
        bamfn: path to BAM file (sorted and indexed)
        fastafn: path to FASTA file (indexed) for the reference genome
        min_mapq: minimum mapping quality (default: 10)
        min_read_quality: minimum mean read quality (default: 5)
        min_base_quality: minimum base quality (default: 5)
        min_depth: minimum read depth (default: 10)
        """
        self.bamfn = bamfn
        self.fastafn = fastafn
        self.seq_params = {
            "min_mapq": min_mapq,
            "min_read_quality": min_read_quality,
            "min_base_quality": min_base_quality,
            "min_depth": min_depth
        }

    def iter(self, chrom, start=0, end=None):
        return ExtractionIterator(self, chrom,start,end)

class ExtractionIterator:

    def __init__(self, extractor, chrom, start=0, end=None):
        assert isinstance(extractor, Extractor)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.bamfile = pysam.AlignmentFile(extractor.bamfn, "rb")
        self.fastafile = pysam.FastaFile(extractor.fastafn)
        self.seq_params = extractor.seq_params
        self.columns = None

    def __iter__(self):
        return self

    def __next__(self):
        if self.columns is None:
            self.columns = self.bamfile.pileup(self.chrom, self.start, self.end, ignore_overlaps=True)
        for column in self.columns:
            site_reads = []
            if column.nsegments < self.seq_params['min_depth']: continue
            pos = column.reference_pos
            if pos < self.start: continue
            if not self.end is None and pos > self.end: break
            if self.fastafile.fetch(self.chrom,pos,pos+2).upper() == "CG":
                site_reads = extract(column,True,self.seq_params)
            elif pos > 0 and self.fastafile.fetch(self.chrom,pos-1,pos+1).upper() == "CG":
                site_reads = extract(column,False,self.seq_params)
            if len(site_reads)> 0: 
                return site_reads
        raise StopIteration


def extract(column,is_fwd,params):
    site_reads = dict()
    for read in column.pileups:
        if read.is_del or read.is_refskip:
            continue
        site_read = SiteRead(read,is_fwd)
        if not site_read.is_good(params): 
            continue
        previous = site_reads.get(site_read.get_name())
        if not previous is None \
           and previous.get_baseq() > site_read.get_baseq(): 
            continue
        site_reads[site_read.get_name()] = site_read
    if len(site_reads) < params['min_depth']:
        return []
    return list(site_reads.values())

