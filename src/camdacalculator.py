import pandas as pd
from collections import deque 

from .methylseqdataset import MethylSeqDataset

class CAMDACalculator:

    def __init__(self,dataset,size):
        assert isinstance(dataset, MethylSeqDataset)
        self.dataset = dataset
        self.window = CAMDAWindow(dataset)
        self.size = size

    def calculate(self,chrom,start=0,end=None):
        for column in self.window.slide(chrom,start,end):
            concurrence_cytosines = 0
            meth = 0
            unmeth = 0
            meth_clones = 0
            unmeth_clones = 0
            concurrence_clones = 0
            region_start = column[0].pos - self.size//2
            region_end = region_start + self.size 
            for column_cread in column:
                clone_meth = 0
                clone_unmeth = 0
                region_meth = 0
                region_unmeth = 0
                for cread in column_cread.read.creads:
                    in_region = cread.pos >= region_start \
                                and cread.pos < region_end
                    if cread.is_methylated:
                        clone_meth += 1
                        if in_region: region_meth += 1
                    else:
                        clone_unmeth += 1
                        if in_region: region_unmeth += 1
                meth += region_meth
                if clone_unmeth > 0:
                    if clone_meth > 0:
                        concurrence_cytosines += region_unmeth
                        concurrence_clones += 1
                    else:
                        unmeth += region_unmeth
                        unmeth_clones += 1
                elif clone_meth > 0:
                    meth_clones += 1
            number_cytosines = concurrence_cytosines + meth + unmeth
            number_clones = concurrence_clones + meth_clones + unmeth_clones
            if number_clones > 0:
                yield { 
                    "chrom": column[0].chrom,
                    "start": region_start,
                    "end": region_end,
                    "depth": number_clones,
                    "nconcurrence": concurrence_cytosines,
                    "camda": concurrence_cytosines/float(number_cytosines),
                    "nconcurrence_clones": concurrence_clones,
                    "nmeth_clones": meth_clones,
                    "unweighted":concurrence_clones/float(number_clones),
                    "nmeth": meth,
                    "meth_pct": meth/float(number_cytosines)
                }


class CAMDAWindow:
    def __init__(self,dataset):
        assert isinstance(dataset,MethylSeqDataset)
        self.dataset = dataset

    def slide(self,chrom,start=0,end=None):
        ## pre_start includes all cytosines covered by reads covering start
        pre_start = start
        iterator = self.dataset.methylation(chrom,start,end)
        for column in iterator:
            for cread in column:
                if cread.read.start < pre_start:
                    pre_start = cread.read.start
            break
        ## load pre_start cytosine data
        iterator = self.dataset.methylation(chrom,pre_start)
        columns = deque()
        for column in iterator:
            if column[0].pos >= start:
                columns.append(column)
                break
        ## begin at first cytosine at or after start
        current_idx = 0
        while current_idx < len(columns):
            if not end is None and columns[current_idx][0].pos > end:
                break
            ## identify boundary of the current window
            window_start = columns[current_idx][0].read.start
            window_end = columns[current_idx][0].read.end
            for cread in columns[current_idx]:
                window_start = min(window_start,cread.read.start)
                window_end = max(window_end,cread.read.end)
            ## cover all cytosines in the window
            for column in iterator:
                columns.append(column)
                if column[0].pos > window_end:
                    break
            ## remove all cytosines outside the window
            while columns[0][0].pos < window_start:
                columns.popleft()
                current_idx -= 1
            yield columns[current_idx]
            current_idx += 1 
