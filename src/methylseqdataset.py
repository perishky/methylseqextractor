import pysam
import numpy
from .cytosineread import CytosineRead
from .read import Read
from .column import Column

class MethylSeqDataset:
    def __init__(
            self, 
            bamfn,
            fastafn,
            merge_strands=True,
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
        self.bamfile = pysam.AlignmentFile(bamfn,"rb")
        self.fastafile = pysam.FastaFile(fastafn)
        self.min_mapq=min_mapq 
        self.min_read_quality=min_read_quality
        self.min_base_quality=min_base_quality
        self.min_depth=min_depth
        self.merge_strands = merge_strands

    def methylation(self,chrom,start=0,end=None):
        if not self.merge_strands:
            return self.methylation_stranded(chrom,start,end)
        else:
            return self.methylation_unstranded(chrom,start,end)

    def methylation_unstranded(self,chrom,start=0,end=None):
        previous_column = None
        for column in self.methylation_stranded(chrom,start,end):
            if column.strand == "-":
                column.pos -= 1
                for cread in column.values():
                    cread.pos -= 1
            if not previous_column is None:
                if previous_column.pos == column.pos:
                    column.merge(previous_column)
                    yield column
                    previous_column = None
                    continue
                else:
                    yield previous_column
                    previous_column = None
            if column.strand == "-":
                yield column
            else:
                previous_column = column
        if not previous_column is None:
            yield previous_column

    def methylation_stranded(self,chrom,start=0,end=None):
        reads = {}
        pileup = self.bamfile.pileup(chrom,start,end,ignore_overlaps=True)
        for pileup_column in pileup:
            if pileup_column.nsegments < self.min_depth: continue
            pos = pileup_column.reference_pos
            if pos < start: continue
            if not end is None and pos > end: break
            if self.fastafile.fetch(chrom,pos,pos+2).upper() == "CG":
                strand = "+"
                site_reads = self.extract_site_reads(pileup_column,True)
            elif pos > 0 and self.fastafile.fetch(chrom,pos-1,pos+1).upper() == "CG":
                strand = "-"
                site_reads = self.extract_site_reads(pileup_column,False)
            else:
                continue
            if len(site_reads)> 0: 
                current_column = Column(chrom,pos,strand)
                for name in site_reads:
                    site_read = site_reads.get(name)
                    base = site_read_base(site_read)
                    baseq = site_read_base_quality(site_read)
                    is_methylated = site_read_is_methylated(site_read)
                    if name in reads:
                        read = reads[name]
                    else:
                        read_start = site_read_start(site_read)
                        read_end = site_read_end(site_read)
                        readq = site_read_quality(site_read)
                        read = Read(name,chrom,read_start,read_end,strand,readq)
                        reads[name] = read
                    cread = CytosineRead(read,current_column,chrom,pos,base,baseq,is_methylated)
                for name in list(reads.keys()):
                    if reads[name].end < pos:
                        del reads[name]
                yield current_column
                
    def extract_site_reads(self,pileup_column,is_fwd):
        site_reads = dict()
        for site_read in pileup_column.pileups:
            if site_read.is_del or site_read.is_refskip:
                continue
            if not self.is_good_site_read(site_read,is_fwd): 
                continue
            previous = site_reads.get(site_read_name(site_read))
            if not previous is None:
                previous_baseq = site_read_base_quality(previous)
                current_baseq = site_read_base_quality(site_read)
                if previous_baseq > current_baseq:
                    continue
            site_reads[site_read.alignment.query_name] = site_read
        if len(site_reads) < self.min_depth:
            return {}
        return site_reads

    def is_good_site_read(self,site_read,is_fwd):
        assert isinstance(site_read, pysam.PileupRead)
        return not (site_read.is_del or site_read.is_refskip) \
            and site_read.alignment.is_proper_pair \
            and site_read_quality(site_read) > self.min_mapq \
            and site_read_quality(site_read) != 255 \
            and site_read_quality(site_read) > self.min_read_quality \
            and site_read_base_quality(site_read) > self.min_base_quality \
            and is_strand_consistent(site_read,is_fwd) \
            and not site_read_methylation_call(site_read) is None

def site_read_name(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    return site_read.alignment.query_name

def site_read_chrom(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    return site_read.alignment.reference_name

def site_read_start(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    return site_read.alignment.reference_start

def site_read_end(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    return site_read.alignment.reference_end

def site_read_pos(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    align = site_read.alignment
    return align.get_reference_positions(full_length=True)[site_read.query_position]

def site_read_base(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    align = site_read.alignment
    return align.query_sequence[site_read.query_position].upper()

def site_read_quality(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    align = site_read.alignment
    return numpy.mean(align.query_qualities)

def site_read_base_quality(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    align = site_read.alignment
    return align.query_qualities[site_read.query_position]

def site_read_methylation_call(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    ## C code from MethylDackel converted to python
    base = site_read.alignment.query_sequence[site_read.query_position].upper()
    cstrand = site_read_conversion_strand(site_read)
    if (cstrand=="OT") or (cstrand=="CTOT"):
        if base=="C": 
            return "m"
        elif base == "T":
            return "u"
    elif (cstrand=="OB") or (cstrand=="CTOB"):
        if base=="G":
            return "m"
        elif base=="A":
            return "u"
    return None

def site_read_conversion_strand(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    ## C code from MethylDackel converted to python
    align = site_read.alignment
    xg = None
    if align.has_tag("XG"):
        xg = align.get_tag("XG")
        if xg != "CT" and xg != "GA": 
            xg = None
    if xg is None: ## BAM file from MethylDackel
        if align.is_proper_pair: 
            if   (align.flag & 0x50) == 0x50: return "OB"  ## read1 -
            elif (align.flag & 0x40): return "OT"          ## read1 +
            elif (align.flag & 0x90) == 0x90: return "OT"  ## read2 -
            elif (align.flag & 0x80): return "OB"          ## read2 +
            return None
        else:
            if (align.flag & 0x10): return "OB"            ## single -
            return "OT"                                    ## single +
    else:  ## BAM file from Bismark
        if xg == "CT":
            if   (align.flag & 0x51) == 0x41: return "OT"  ## read1 +
            elif (align.flag & 0x51) == 0x51: return "CTOT"## read1 -
            elif (align.flag & 0x91) == 0x81: return "CTOT"## read2 +
            elif (align.flag & 0x91) == 0x91: return "OT"  ## read2 -
            elif (align.flag & 0x10): return "CTOT"        ## single -
            return "OT"                                    ## single + 
        else:
            if (align.flag & 0x51) == 0x41: return "CTOB"  ## read1 + 
            elif (align.flag & 0x51) == 0x51: return "OB"  ## read1 - 
            elif (align.flag & 0x91) == 0x81: return "OB"  ## read2 + 
            elif (align.flag & 0x91) == 0x91: return "CTOB"## read2 - 
            elif (align.flag & 0x10): return "OB"          ## single - 
            return "CTOB"                                  ## single +

def site_read_is_methylated(site_read):
    assert isinstance(site_read, pysam.PileupRead)
    return site_read_methylation_call(site_read) == "m"

def is_strand_consistent(site_read,is_fwd):
    assert isinstance(site_read, pysam.PileupRead)
    cstrand = site_read_conversion_strand(site_read)
    return is_fwd and (cstrand in ['OT','CTOT']) \
        or not is_fwd and (cstrand in ['OB','CTOB'])

