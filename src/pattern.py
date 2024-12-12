import pandas as pd

class Pattern:
    def __init__(self,chrom,positions,reads,meth):
        assert isinstance(reads,list)
        assert isinstance(meth,list)
        assert len(reads) == len(meth)
        assert len(positions) == len(meth[0])
        self.chrom=chrom
        self.positions=positions
        self.reads=reads
        self.meth=meth

    def __str__(self):
        meth = self.meth
        meth = [["" if m < 0 else str(m) for m in clone] for clone in meth]
        df = pd.DataFrame(meth)
        df.insert(0,"read",self.reads)
        out = ("chrom = " + self.chrom 
                + " positions = " + str(self.positions) + "\n" 
                +  df.to_string(header=False,index=False))
        return out
