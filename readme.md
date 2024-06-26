# methylseqextractor

A tool allowing the flexible calculation of concordance
and heterogeneity metrics in methyl-seq data.

## Installing dependencies

```
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

```
cd genome
samtools faidx hg19.fa
```

Index the input BAM file.

```
cd data
samtools index sample.bam
```

## Extracting reads for all CpG sites

Prepare the environment


Run the analysis
```
mkdir output-reads
python extract.py \
  data/sample.bam \
  genome/hg19.fa \
  output-reads/methylation_per_read
```

Output is a csv file for each chromosome (`output/methylation_per_read_chr*.csv`)
listing the reads that overlap each CpG site (with at least 10 reads).

With 12 processors, should take about 1-2 minutes. 


## Calculating methylation levels

```
mkdir output-levels
python calculate-levels.py \
  dat/sample.bam \
  genome/hg19.fa \
  output-levels/methylation_levels
```

## Calculating heterogeneity


