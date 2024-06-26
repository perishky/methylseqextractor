import pysam
import numpy

class MethylSeqExtractor:

    """Iterates through methyl-seq reads by CpG site"""

    def __init__(
        self, 
        bamfn,
        fastafn,
        chrom,
        start=0,
        stop=None,
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
        chrom: chromosome
        start: chromosomal start position (default: 0)
        stop: chromosomal end position (default: None)
        min_mapq: minimum mapping quality (default: 10)
        min_read_quality: minimum mean read quality (default: 5)
        min_base_quality: minimum base quality (default: 5)
        min_depth: minimum read depth (default: 10)

        Returns 
        -------
        List of dictionaries, one per read, intersecting the current 
        CpG site including: 
        - read query name ('name'), 
        - read strand of converted genome ('strand'), 
        - read mapping quality ('mapq'), and
        - mean quality of bases in the read ('readq'), 
        - chromosome ('chrom'), 
        - chromosomal position of the CpG site ('pos', 1-based), 
        - sequenced base ('base'), 
        - quality of the base('baseq') and 
        - whether not it was called methylated ('is_methylated').
        """
        self.bamfile = pysam.AlignmentFile(bamfn, "rb")
        self.fastafile = pysam.FastaFile(fastafn)
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.seq_params = {
            "min_mapq": min_mapq,
            "min_read_quality": min_read_quality,
            "min_base_quality": min_base_quality,
            "min_depth": min_depth
        }
        self.columns = None

    def __iter__(self):
        return self

    def __next__(self):
        if self.columns is None:
            self.columns = self.bamfile.pileup(self.chrom, self.start, self.stop, ignore_overlaps=True)
        for column in self.columns:
            sites = []
            if column.nsegments < self.seq_params['min_depth']: continue
            pos = column.reference_pos
            if self.fastafile.fetch(self.chrom,pos,pos+2).upper() == "CG":
                sites = extract_cpg_site_methylation(column,True,self.seq_params)
            elif pos > 0 and self.fastafile.fetch(self.chrom,pos-1,pos+1).upper() == "CG":
                sites = extract_cpg_site_methylation(column,False,self.seq_params)
            if len(sites)==0: continue
            return sites
        raise StopIteration

def is_good_read(read,params):
    align = read.alignment
    return not (read.is_del or read.is_refskip) \
        and align.is_proper_pair \
        and align.mapping_quality > params['min_mapq'] \
        and align.mapping_quality != 255 \
        and numpy.mean(align.query_qualities) > params['min_read_quality'] \
        and align.query_qualities[read.query_position] > params['min_base_quality']

def extract_cpg_site_methylation(column,is_fwd,params):
    sites = dict()
    pos = column.reference_pos
    for read in column.pileups:
        align = read.alignment
        if not is_good_read(read,params): continue
        baseq = align.query_qualities[read.query_position]
        if (align.query_name in sites) and baseq > sites[align.query_name]['baseq']: continue
        strand = get_conversion_strand(read)
        if is_fwd and strand != "OT" and strand != "CTOT": continue
        if not is_fwd and strand != "OB" and strand != "CTOB": continue
        call = methylation_call(read)
        if call is None: continue
        base = align.query_sequence[read.query_position].upper()
        sites[align.query_name] = { 
            "name": align.query_name,
            "first": align.is_read1,
            "chrom": column.reference_name,
            "pos": pos+1,
            "strand": strand,
            "base": base,
            "depth": column.nsegments,
            "mapq": align.mapping_quality,
            "readq": numpy.mean(align.query_qualities),
            "baseq": baseq,
            "is_methylated": call=="m"
        }
    ## skip CpG site if depth is too low
    if len(sites) < params['min_depth']:
        return []
    return sites.values()

def methylation_call(read):
    ## C code from MethylDackel converted to python
    align = read.alignment
    base = align.query_sequence[read.query_position].upper()
    strand = get_conversion_strand(read)
    if (strand=="OT") or (strand=="CTOT"):
        if base=="C": 
            return "m"
        elif base == "T":
            return "u"
    elif (strand=="OB") or (strand=="CTOB"):
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
        
        
