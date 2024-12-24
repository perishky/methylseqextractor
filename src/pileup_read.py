import pysam
import numpy

def is_good(pileup_read,strand,min_mapq,min_read_quality,min_base_quality):
    assert isinstance(pileup_read, pysam.PileupRead)        
    return not (pileup_read.is_del or pileup_read.is_refskip) \
        and pileup_read.alignment.is_proper_pair \
        and get_quality(pileup_read) > min_mapq \
        and get_quality(pileup_read) != 255 \
        and get_quality(pileup_read) > min_read_quality \
        and get_base_quality(pileup_read) > min_base_quality \
        and is_strand_consistent(pileup_read,strand) \
        and not get_methylation_call(pileup_read) is None

def get_read_name(pileup_read,strand):
    assert isinstance(pileup_read, pysam.PileupRead)
    align = pileup_read.alignment
    return (
        align.query_name 
        + "-" + str(align.reference_name)
        + ":" + str(align.reference_start)
        + "-" + str(align.reference_end)
        + strand)

def get_chrom(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    return pileup_read.alignment.reference_name

def get_start(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    return pileup_read.alignment.reference_start

def get_end(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    return pileup_read.alignment.reference_end

def get_pos(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    align = pileup_read.alignment
    return align.get_reference_positions(full_length=True)[pileup_read.query_position]

def get_base(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    align = pileup_read.alignment
    return align.query_sequence[pileup_read.query_position].upper()

def get_quality(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    align = pileup_read.alignment
    return numpy.mean(align.query_qualities)

def get_base_quality(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    align = pileup_read.alignment
    return align.query_qualities[pileup_read.query_position]

def get_methylation_call(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    ## C code from MethylDackel converted to python
    base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
    cstrand = get_conversion_strand(pileup_read)
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

def get_conversion_strand(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    ## C code from MethylDackel converted to python
    align = pileup_read.alignment
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

def is_methylated(pileup_read):
    assert isinstance(pileup_read, pysam.PileupRead)
    return get_methylation_call(pileup_read) == "m"

def is_strand_consistent(pileup_read,strand):
    assert isinstance(pileup_read, pysam.PileupRead)
    cstrand = get_conversion_strand(pileup_read)
    return ((strand == "+") and (cstrand in ['OT','CTOT'])
            or (strand == "-") and (cstrand in ['OB','CTOB']))

