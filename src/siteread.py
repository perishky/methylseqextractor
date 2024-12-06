import pysam
import numpy

class SiteRead (dict):
    """
        - name: read query name  
        - first: is the first in a pair of reads 
        - chrom: chromosome
        - pos: chromosomal position of the CpG site (0-based)
        - base: sequenced base 
        - strand: forward or reverse strand
        - cstrand: alignment strand of converted genome
        - mapq: read mapping quality 
        - readq: mean quality of bases in the read
        - baseq: quality of the base
        - is_methylated: whether not it was called methylated 
    """

    def __init__(self,read,is_fwd): 
        assert isinstance(read, pysam.PileupRead)
        self.read = read
        align = read.alignment
        super().__init__(
            name=align.query_name,
            first=align.is_read1,
            chrom=align.reference_name,
            pos=align.get_reference_positions(full_length=True)[read.query_position],
            strand=("+" if is_fwd else "-"),
            cstrand=get_conversion_strand(read),
            base=align.query_sequence[read.query_position].upper(),
            mapq=align.mapping_quality,
            readq=numpy.mean(align.query_qualities),
            baseq=align.query_qualities[read.query_position],
            is_methylated=(methylation_call(read)=="m"))

    def get_name(self): return self['name']
    def get_chrom(self): return self['chrom']
    def get_strand(self): return self['strand']
    def get_pos(self,merge_strands=False):
        if self['strand'] == "+" or not merge_strands:
            return self['pos']
        else:
            return self['pos'] - 1
    def get_baseq(self): return self['baseq']
    def is_methylated(self): return self['is_methylated']
    def is_strand_consistent(self):
        return ((self['strand'] == "+") \
                and (self['cstrand'] in['OT','CTOT'])) \
            or ((self['strand'] == "-") \
                and (self['cstrand'] in['OB','CTOB']))

    def is_good(self,params):
        return not (self.read.is_del or self.read.is_refskip) \
            and self.read.alignment.is_proper_pair \
            and self['mapq'] > params['min_mapq'] \
            and self['mapq'] != 255 \
            and self['readq'] > params['min_read_quality'] \
            and self['baseq'] > params['min_base_quality'] \
            and self.is_strand_consistent() \
            and not self.is_methylated() is None

def methylation_call(read):
    ## C code from MethylDackel converted to python
    align = read.alignment
    base = align.query_sequence[read.query_position].upper()
    cstrand = get_conversion_strand(read)
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

def get_conversion_strand(read):
    ## C code from MethylDackel converted to python
    align = read.alignment
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
        
