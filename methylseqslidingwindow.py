from methylseqextractor import MethylseqExtractor
from collections import deque 
import numpy


class MethylSeqSlidingWindow: 

    """Slides a window across the genome showing methylation patterns"""

    def __init__(self, window, extractor):
        """
        Parameter
        ---------
        window: integer indicating window size in base-pairs
        extractor: MethylSeqExtractor

        Returns
        -------
        Dictionary for the current window including:
        - chromosome ('chrom')
        - chromosomal positions of CpG sites ('pos':array)
        - methylation states of CpG site sites in each read ('meth')
          (each read is a key in a dictionary)
        """
        self.window = window
        self.extractor = extractor
        self.cytosines = deque()

    def __iter__(self):
        return self:

    def __next__(self):
        for cytosines in self.extractor: 
            if len(self.cytosines) == 0 \
               or in_same_window(cytosines[0],self.cytosines[0],self.window):
                self.cytosines.extend(cytosines)
            else:
                break
        pattern = methylation_pattern(self.cytosines)
        while len(self.cytosines) > 0 \
              and not in_same_window(cytosines[0],self.cytosines[0],self.window):
            self.cytosines.popleft()
        self.cytosines.extend(cytosines)
        return pattern

def in_same_window(cytosine1,cytosine2,window):
    return (cytosine1['chrom'] == cytosine2['chrom']) \
        and (abs(cytosine1['pos']-cytosine2['pos']) < window) 

def methylation_pattern(cytosines):
    reads = dict()
    positions = []
    for cytosine in cytosines:
        reads[cytosine['name']] = None
        if positions[-1] != cytosine['pos']:
            positions.append(cytosine['pos'])
    for name in reads.keys():
        reads[name] = [-1]*len(positions)
    idx = 0
    for cytosine in cytosines:
        while positions[idx] < cytosine['pos']: 
            idx += 1
        reads[cytosine['name']][idx] = 1 if cytosine['is_methylated'] else 0
    {
        chrom: cytosine['chrom'],
        pos: positions,
        meth: reads
    }

