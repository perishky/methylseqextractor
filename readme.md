# methylseqextractor

A tool allowing the flexible calculation of
DNA methylation concordance
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

Each analysis is applied to a single sorted and indexed BAM file.
If necessary, a BAM file can be sorted and indexed using samtools, e.g.

```bash
samtools sort sample_unsorted.bam -o sample.bam
samtools index sample.bam
```

All analyses require the reference genome sequence
in an indexed FASTA file.

If necessary, a FASTA file can be indexed using samtools.

```bash
samtools faidx hg19.fa
```

## Running analyses

See [example.html](example.html) for an example
as well as instructions for how to make use of multiple processors
for efficient genome-wide analyses. 

## Generating example.html

Generating example output file `example.html` from `example.qmd`
uses the `jupyter` Python package.
It can be installed as follows.

```bash
pip3 install jupyter
```

The code assumes the existence of `data/sample.bam` containing sequencing reads. 
A small example file was used to the generate `example.html`.
It was obtained from a real bam file (called `original.bam`)
but down-sampled to simplify the output.

```bash
mkdir data
cd data
## if necessary, install samtools
conda install bioconda::samtools
## select reads from just chromosome 1
samtools view -bo chr1.bam original.bam chr1
## randomly retain 10% of the reads
samtools view -b -s 0.1 chr1.bam > sample.bam
## index the resulting file
samtools index sample.bam
```

The output file `example.html` can be generated using 'quarto'. 

```bash
quarto render example.qmd
```

