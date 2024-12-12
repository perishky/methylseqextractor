## exec(open("example.py").read())

# from methylseqextractor import (
#     Extractor, 
#     WindowMaker, 
#     LevelCalculator, 
#     WindowSlider, 
#     ClonalFlipCounter,
#     CAMDACalculator,
#     ConcurrenceCalculator
# )

from src import (
    Extractor, 
    WindowMaker, 
    LevelCalculator, 
    WindowSlider, 
    ClonalFlipCounter,
    CAMDACalculator,
    ConcurrenceCalculator
)

import pandas

## Extract reads containing specified CpG sites

bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
extractor = Extractor(bamfn, fastafn)
iter = extractor.iter("chr1",1245000,1246000)

site_reads = pandas.DataFrame(next(iter))
site_reads

## Calculate DNA methylation levels

levelcalculator = LevelCalculator(extractor)
iter = levelcalculator.iter("chr1", 1245000, 1246000)
levels = pandas.DataFrame([site for site in iter])

levels[0:5]

## View DNA methylation patterns

slider = WindowSlider(WindowMaker(1000), extractor)
iter = slider.iter("chr1", 1245000, 1246000)
print(next(iter)) ## first window position
 
## Count clonal flips

counter = ClonalFlipCounter(slider)
iter = counter.iter("chr1", 1245000, 1246000)
flips = pandas.DataFrame([ region for region in iter])

flips[0:5]

## Calculate CAMDA scores 

camdacalculator = CAMDACalculator(500,250,extractor)
iter = camdacalculator.iter("chr1", 1245000, 1246000)
scores = pandas.DataFrame([site for site in iter])

scores[0:5]

## Calculate concurrence scores

concurrencecalculator = ConcurrenceCalculator(slider)
iter = concurrencecalculator.iter("chr1", 1245000, 1246000)
scores = pandas.DataFrame([site for site in iter])

scores[0:5]

## Analysing the entire genome efficiently using multiple processors

import multiprocessing

def calculate_chrom_concurrences(chrom):
  regions = concurrencecalculator.iter(chrom)
  return pandas.DataFrame([region for region in regions])

chromosomes = ["chr"+str(i) for i in range(1,23)] + ["chrX"]

with multiprocessing.Pool(processes=12) as pool:
    for stats in pool.imap(calculate_chrom_concurrences, chromosomes):
        if len(stats) > 0:
            stats.to_csv("clonalflipcounts_" + stats['chrom'][0] + ".csv")
