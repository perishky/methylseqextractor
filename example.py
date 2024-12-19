## exec(open("example.py").read())

from src import (
  Window,
  MethylSeqDataset,
  LevelCalculator,
  ClonalFlipCounter,
  ConcurrenceCalculator,
  CAMDACalculator
)

import pandas

## Extract reads containing specified CpG sites

bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
dataset = MethylSeqDataset(bamfn, fastafn)
iter = dataset.methylation("chr1",1245000,1246000)

sites = next(iter)
pandas.DataFrame([site.get_dict() for site in sites])

## Calculate DNA methylation levels

levels = LevelCalculator(dataset)
iter = levels.calculate("chr1", 1245000, 1246000)
scores = pandas.DataFrame([site for site in iter])

scores[0:5]

## View DNA methylation patterns

window = Window(dataset,200)
iter = window.slide("chr1", 1245000, 1246000)
print(next(iter)) ## first window position
print(next(iter)) ## second window position
print(next(iter)) ## third window position
 
## Count clonal flips

flips = ClonalFlipCounter(dataset,200)
iter = flips.calculate("chr1", 1245000, 1246000)
scores = pandas.DataFrame([ region for region in iter])

scores[0:5]

## Calculate CAMDA scores 
camda = CAMDACalculator(dataset, 1)
iter = camda.calculate("chr1", 1245000, 1246000)
scores = pandas.DataFrame([site for site in iter])

scores[0:5]

## Calculate concurrence scores

concurrence = ConcurrenceCalculator(dataset,200)
iter = concurrence.calculate("chr1", 1245000, 1246000)
scores = pandas.DataFrame([region for region in iter])

scores[0:5]

## Analysing the entire genome efficiently using multiple processors

if False:
  import multiprocessing
  concurrence = ConcurrenceCalculator(dataset,1000)

  def calculate_chrom_concurrences(chrom):
    iter = concurrence.calculate(chrom)
    return pandas.DataFrame([region for region in iter])

  chromosomes = ["chr"+str(i) for i in range(1,23)] + ["chrX"]

  with multiprocessing.Pool(processes=12) as pool:
    for stats in pool.imap(calculate_chrom_concurrences, chromosomes):
      if len(stats) > 0:
        stats.to_csv("concurrences_" + stats['chrom'][0] + ".csv")
