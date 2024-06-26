# methylseqextractor

A tool allowing the flexible calculation of concordance
and heterogeneity metrics in methyl-seq data.

## Installing dependencies

```{bash}
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

```{bash}
cd genome
samtools faidx hg19.fa
```

Index the input BAM file.

```{bash}
cd data
samtools index sample.bam
```

## Extract reads containing specified CpG sites

Apply to just CpG sites in the region chr1:1245000-1246000

```{python}
from methylseqextractor import MethylSeqExtractor
import pandas
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
extractor = MethylSeqExtractor(bamfn, fastafn, "chr1", 1245000, 1246000)
pandas.DataFrame([site_read for site_read in extractor])
```

Apply to the entire genome

```{bash}
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

```{python}
from methylseqextractor import MethylSeqExtractor
from methylseqlevels import MethylSeqLevels
import pandas
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
extractor = MethylSeqExtractor(bamfn, fastafn, "chr1", 1245000, 1246000)
pandas.DataFrame([site for site in MethylSeqLevels(extractor)])
```

Apply across the entire genome

```{bash}
mkdir output-levels
python calculate-levels.py \
  dat/sample.bam \
  genome/hg19.fa \
  output-levels/methylation_levels
```

## Slide DNA methylation patterns window across the genome

```{python}
from methylseqextractor import MethylSeqExtractor
from methylseqslidingwindow import MethylSeqSlidingWindow
import pandas
bamfn = "data/sample.bam"
fastafn = "genome/hg19.fa"
window = 200 ## base pairs
extractor = MethylSeqExtractor(bamfn, fastafn, "chr1", 1245000, 1246000)

for window in MethylSeqSlidingWindow(window,extractor):
    start = window['pos'][0]
    end = window['pos'][-1]
    print(window['chrom'] + ":", str(start), "-", str(end))
    for read in window['meth']:
        print(read)
```


## Calculate CAMDA scores across the genome
