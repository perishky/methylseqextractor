import pysam
from .cytosineread import CytosineRead
from .read import Read
from .column import Column
from . import pileup_read

class MethylSeqDataset:
    def __init__(
            self, 
            bamfn,
            fastafn,
            merge_strands=True,
            min_mapq=10, 
            min_read_quality=5,
            min_base_quality=5
    ):
        """
        Parameters
        ----------
        bamfn: path to BAM file (sorted and indexed)
        fastafn: path to FASTA file (indexed) for the reference genome
        min_mapq: minimum mapping quality (default: 10)
        min_read_quality: minimum mean read quality (default: 5)
        min_base_quality: minimum base quality (default: 5)
        """
        self.bamfile = pysam.AlignmentFile(bamfn,"rb")
        self.fastafile = pysam.FastaFile(fastafn)
        self.min_mapq=min_mapq 
        self.min_read_quality=min_read_quality
        self.min_base_quality=min_base_quality
        self.merge_strands = merge_strands

    def methylation(self,chrom,start=0,end=None):
        if self.merge_strands:
            return self._methylation_unstranded(chrom,start,end)
        else:
            return self._methylation_stranded(chrom,start,end)

    def _methylation_unstranded(self,chrom,start=0,end=None):
        previous_column = None
        for column in self._methylation_stranded(chrom,start,end):
            if not previous_column is None:
                if previous_column.get_pos() == column.get_pos():
                    previous_column.merge_negative(column)
                    yield previous_column
                    previous_column = None
                    continue
                else:
                    yield previous_column
                    previous_column = None
            previous_column = column
        if not previous_column is None and previous_column.get_pos() < end:
            yield previous_column

    def _methylation_stranded(self,chrom,start=0,end=None):
        reads = dict()
        pileup = self.bamfile.pileup(chrom,start,end,ignore_overlaps=True)
        for pysam_column in pileup:
            pos = pysam_column.reference_pos
            if pos < start: continue
            if not end is None and pos > end-1: break
            if self.fastafile.fetch(chrom,pos,pos+2).upper() == "CG":
                strand = "+"
            elif pos > 0 and self.fastafile.fetch(chrom,pos-1,pos+1).upper() == "CG":
                strand = "-"
            else:
                continue
            column = Column(chrom,pos,strand,self.merge_strands)
            for pread in pysam_column.pileups:
                is_good = pileup_read.is_good(
                    pread,
                    strand,
                    self.min_mapq,
                    self.min_read_quality,
                    self.min_base_quality)
                if not is_good:
                    continue
                name = pileup_read.get_read_name(pread,strand)
                if name in reads:
                    read = reads[name]
                else:
                    read = Read(
                        name,
                        chrom,
                        pileup_read.get_start(pread),
                        pileup_read.get_end(pread),
                        strand,
                        pileup_read.get_quality(pread),
                        self.merge_strands)
                    reads[name] = read
                cread = CytosineRead(
                    read,
                    pos,
                    pileup_read.get_base(pread),
                    pileup_read.get_base_quality(pread),
                    pileup_read.is_methylated(pread),
                    self.merge_strands)
                read.add_cread(cread) 
                column.add_cread(cread)
            for name in list(reads.keys()):
                if reads[name].get_end() < pos+1:
                    del reads[name]
            yield column
                
