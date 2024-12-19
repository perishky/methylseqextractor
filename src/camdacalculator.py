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
        column = next(self.dataset.methylation(chrom,start,end))
        virtual_start = start
        for cread in column:
            if cread.read.start < virtual_start:
                virtual_start = cread.read.start
        for column in self.window.slide(chrom,virtual_start,end):
            if column[0].pos < start:
                continue
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
                    in_region = cread.pos >= region_start and cread.pos < region_end
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
                    "nconcurrence": concurrence_cytosines,
                    "camda": concurrence_cytosines/float(number_cytosines),
                    "nconcurrence_clones": concurrence_clones,
                    "nmeth_clones": meth_clones,
                    "unweighted_camda":concurrence_clones/float(number_clones),
                    "nmeth": meth,
                    "meth_pct": meth/float(number_cytosines)
                }


class CAMDAWindow:
    def __init__(self,dataset):
        assert isinstance(dataset,MethylSeqDataset)
        self.dataset = dataset

    def slide(self,chrom,start=0,end=None):
        iterator = self.dataset.methylation(chrom,start,end)
        columns = deque()
        current_idx = -1
        next_column = None
        try:
            columns.append(next(iterator))
            next_column = next(iterator)
        except StopIteration:
            if len(columns) == 0:
                return
        for column in iterator:
            current_idx += 1 ## current position = columns[current_idx][0].pos
            start = min([cread.read.start for cread in columns[current_idx]])
            end = max([cread.read.end for cread in columns[current_idx]])
            ## add columns with cytosines in reads covering the current position
            while not next_column is None and next_column[0].pos <= end:
                columns.append(next_column)
                try:
                    next_column = next(iterator)
                except StopIteration:
                    next_column = None
                    break
            ## remove columns for cytosines not in any reads covering the current position
            while columns[0][0].pos < start:
                columns.popleft()
                current_idx -= 1
                if len(columns) == 0:
                    return
            yield columns[current_idx]
