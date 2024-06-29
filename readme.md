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

## Extract reads containing specified CpG sites

Apply to just CpG sites in the region chr1:1245000-1246000

```python
from extractor import Extractor
import pandas
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
extractor = Extractor(bamfn, fastafn)
pandas.DataFrame([site_read for site_read in extractor.iter("chr1",1245000,1246000])
```

Apply to the entire genome

```bash
mkdir output-reads
python extract.py \
  data/sample.bam \
  genome/hg19.fa \
  output-reads/methylation_per_read
```

Output will appear in a csv file for each chromosome
(`output/methylation_per_read_chr*.csv`)
listing the reads that overlap each CpG site (with at least 10 reads).
With 12 processors, should take about 1-2 minutes. 


## Calculate methylation levels

Calculate methylation levels in just the region chr1:1245000-1246000

```python
from extractor import Extractor
from levelcalculator import LevelCalculator
import pandas
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
calculator = LevelCalculator(Extractor(bamfn, fastafn))
pandas.DataFrame([site for site in calculator.iter("chr1", 1245000, 1246000)])
```

Apply across the entire genome

```bash
mkdir output-levels
python calculate-levels.py \
  dat/sample.bam \
  genome/hg19.fa \
  output-levels/methylation_levels
```

## Slide DNA methylation patterns window across the genome

```python
from extractor import Extractor
from windowslider import WindowSlider
from windowmaker import WindowMaker
import pandas
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
size = 200 ## base pairs
slider = WindowSlider(WindowMaker(size), Extractor(bamfn, fastafn))

for window in slider.iter("chr1", 1245000, 1246000):
    start = window['positions'][0]
    end = window['positions'][-1]
    print(window['chrom'] + ":", str(start), "-", str(end))
    for read in window['meth']:
        print(read)
```


## Calculate clonal flip scores

```
from extractor import Extractor
from windowslider import WindowSlider
from windowmaker import WindowMaker
import pandas
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
size = 200 ## base pairs
counter = ClonalFlipCounter(WindowSlider(WindowMaker(size), Extractor(bamfn, fastafn)))

pandas.DataFrame([site for site in counter.iter("chr1", 1245000, 1246000)])
```
