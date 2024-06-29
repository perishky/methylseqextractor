# methylseqextractor

A tool allowing the flexible calculation of concordance
and heterogeneity metrics in methyl-seq data.

## Installing dependencies

```bash
conda config --add channels r
conda config --add channels bioconda
conda create -n methex python=3.8
conda activate methex
conda install -c bioconda pysam
pip3 install pandas
```

## Preparing the data


Obtain the reference genome FASTA file
and calculate an index.

```bash
cd genome
samtools faidx hg19.fa
```

Index the input BAM file.

```bash
cd data
samtools index sample.bam
```

## Analysis in python

### Import methylseqextractor modules

```python
from methylseqextractor import (
    Extractor, 
    WindowMaker, 
    LevelCalculator, 
    WindowSlider, 
    ClonalFlipCounter,
    ConcurrenceCalculator
)
```

We'll also use `pandas` for data frames.
```python
import pandas
```

### Extract reads containing specified CpG sites

Extract reads overlapping CpG sites in the region chr1:1245000-1246000

```python
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
extractor = Extractor(bamfn, fastafn)
pandas.DataFrame([site_read for site_read in extractor.iter("chr1",1245000,1246000])
```

### Calculate DNA methylation levels

Calculate DNA methylation levels in just the region chr1:1245000-1246000

```python
levelcalculator = LevelCalculator(extractor)
pandas.DataFrame([site for site in levelcalculator.iter("chr1", 1245000, 1246000)])
```

## Slide window across the genome to view DNA methylation patterns

```python
slider = WindowSlider(WindowMaker(50))
for region in slider.iter("chr1", 1245000, 1246000):
  print(region)
```

## Count clonal flips

```python
counter = ClonalFlipCounter(slider)
pandas.DataFrame([site for site in counter.iter("chr1", 1245000, 1246000)])
```

## Calculate concurrence scores

```python
concurrencecalculator = ConcurrenceCalculator(slider)
pandas.DataFrame([site for site in concurrencecalculator.iter("chr1", 1245000, 1246000)])
```

## Analysing the entire genome

Outputs are likely to get large and computationally more costly
when generated for the entire genome.
Forutnately, these analyses are naturally parallalizable
by chromosome using the `multiprocessing`library.

```python
import multiprocessing

def calculate_chrom_concurrences(chrom):
  regions = concurrencecalculator.iter(chrom)
  return pd.DataFrame([region for region in regions])

chromosomes = ["chr"+str(i) for i in range(1,23)] + ["chrX"]

with multiprocessing.Pool(processes=12) as pool:
  for stats in pool.imap(calculate_chrom_concurrences, chromosomes):
    dat.to_csv("clonalflipcounts_" + stats['chrom'][0] + ".csv")
```
