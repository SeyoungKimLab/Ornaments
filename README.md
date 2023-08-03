# ornaments
__ornaments__ is a wrapper program around kallisto to facilitate
variant-aware pseudoalignment of RNA-Seq data to quantify allele-specific
expression or obtain allele-specific read counts. __ornaments__ takes as input
a list of variants (both SNPs and indels are supported) and a reference
transcriptome to create a modified reference that can accommodate allele-specific
alignments across different samples. __ornaments__ requires that the variants
are listed in terms of transcriptomic coordinates. We provide helper scripts
to prepare the ornaments reference as well as transform the variant coordinates.

## Setup
To compile ornaments, cmake is required. First, clone the repository and go into the
ornamnets folder. Once in the ornaments folder, execute the following commands.

```
mkdir build
cd build; cmake ..; make; cd ..
```

An executable should be built with the name _build/src/ornaments_.

## Data
Example ranscriptome files and gene annotation files necessary for the scripts can be found at: [here](https://www.gencodegenes.org/human/)

You can download a list of variant files [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/). The files of interest are the .vcf.gz files which contain a list of all the phased variants, separated by chromosome.

## Workflow
To begin, we must translate our set of genomic variants into transcriptomic variants. 
The provided script `convert_genomic_vcf_to_transcriptomic_vcf.py` will transform
genomic variants (one chromosome at a time, in vcf format) to transcriptomic variants.

Note that you must have [biopython](https://anaconda.org/conda-forge/biopython) and [gtfparse](https://anaconda.org/bioconda/gtfparse) installed in your conda environment.

Example Usage:

```
python3 convert_genomic_vcf_to_transcriptomic_vcf.py -v all.chr1.vcf.gz -g grch38.annotation.gtf -t gencode.v36.transcripts.fa -o chr1.transcriptome.vcf -c chr1.transcriptome.coords
```


Here, `all.chr1.vcf.gz`is the genome VCF file, `grc38.annotation.gtf` is the gene to
transcript annotation file, `gencode.v36.transcripts.fa` is the transcriptome, 
`chr1.transcriptome.vcf` is the output file for the converted variants, and `chr1.transcriptome.coords` 
is an output file that keeps track of the conversions. Please note that a compatible
reference transcriptome and transcript annotation file should be used. If you are getting both files from the same source, make sure that the version is the same.

Upon creating the transcriptomic variants, we need to generate a modified reference
transcriptome to use with the ornaments software. To do this, we need to merge the files
the transcriptome variant files for each chromosome generated in the previous step. The 
`merge_transcriptome_vcfs.py` script is provided for this.

Example Usage:
```
python3 merge_transcriptome_vcfs.py -i in_dir -o transcriptome.vcf -s sample_file
```

Here, `in_dir` is the input directory, which should contain the transcriptome vcf files 
(which are assumed to be named in the format `chr*.transcriptome.vcf`), `transcriptome.vcf` is 
the output file, and `sample_file` is a file with a list of sample names (separated by newlines)
which are considered, since only variants that are heterozygous in at least one of the samples 
will be retained. Alternatively, one can specify a single sample in place of a sample file.

Upon merging, we can use the provided script `create_personalized_transcriptome.py` to modify the reference.

Example Usage to Create Diploid Transcriptome:
```
python3 create_personalized_transcriptome.py -f gencode.v36.transcripts.fa -v chr1.transcriptome.vcf -t diploid -o SAMPLE.diploid -s SAMPLE
```

Here, `diploid` denotes that a diploid transcriptome will be created, `SAMPLE.diploid` denotes the desired filename prefix for the output,
and `SAMPLE` is the name of the sample. The name that should be used is the name that refers to the sample in the VCF file.

Example Usage to Create Ornament Transcriptome:
```
python3 create_personalized_transcriptome.py -f gencode.v36.transcripts.fa -v chr1.transcriptome.vcf -t ornament -o SAMPLE.ornament -s SAMPLE
```

This will result in a file called `SAMPLE.ornament.fa`, which can be used by the ornaments. Ornaments
expands the functionality of __kallisto__, meaning it functions differently only if `SAMPLE.ornament.fa`
is provided as input. Otherwise, it should behave the exact same as __kallisto__.

Running Ornaments:

```
build/src/ornaments index -i SAMPLE.ornament.fa SAMPLE.ornament.index
```

This will construct an ornament index file, called `SAMPLE.ornament.index`, which will be used for pseudoalignment in the quantification. To perform the quantification, we need to sort the `transcriptome.vcf` file that is produced by the merge script. This will be updated in future versions to avoid this explicit sorting step. To run the sorting, you can use the command:

```
cat transcriptome.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > sorted.transcriptome.vcf
```

To perform variant-aware quantification, we need to provide a few more parameters. To run with paired-end data that is specific to a given `SAMPLE`:

```
build/src/ornaments quant -i SAMPLE.ornament.index -o output_dir --vcf sorted.transcriptome.vcf --sample SAMPLE --phased SAMPLE_1.fastq SAMPLE_2.fastq
```

Here, the `vcf` parameter is used to read the variant information, and the `phased` parameter tells ornaments that the data is assumed to be phased. To run ornaments in unphased mode to obtain expected read counts, remove the `phased` parameter. The output of this run mode will be a read counts file.

To do the same with single-end RNA-seq data, we additionally need to specify the fragment length mean and standard deviation used during sequencing:

```
build/src/ quant -i SAMPLE.ornament.index -o output_dir --vcf sorted.transcriptome.vcf --sample SAMPLE --phased --single -l 160 -s 25 SAMPLE_1.fastq SAMPLE_2.fastq
```

In this example the fragment length mean is 160 bp and the standard deviation is 25 bp (you can use these values for Geuvadis data).

Ornaments is built as an addition to kallisto.

# kallisto

__kallisto__ is a program for quantifying abundances of transcripts from
RNA-Seq data, or more generally of target sequences using high-throughput
sequencing reads. It is based on the novel idea of _pseudoalignment_ for
rapidly determining the compatibility of reads with targets, without the need
for alignment. On benchmarks with standard RNA-Seq data, __kallisto__ can
quantify 30 million human bulk RNA-seq reads in less than 3  minutes on a Mac desktop
computer using only the read sequences and a transcriptome index that
itself takes than 10 minutes to build. Pseudoalignment of reads
preserves the key information needed for quantification, and __kallisto__
is therefore not only fast, but also comparably accurate to other existing
quantification tools. In fact, because the pseudoalignment procedure is
robust to errors in the reads, in many benchmarks __kallisto__
significantly outperforms existing tools. The __kallisto__ algorithms are described in more detail in:

NL Bray, H Pimentel, P Melsted and L Pachter, [Near optimal probabilistic RNA-seq quantification](http://www.nature.com/nbt/journal/v34/n5/abs/nbt.3519.html), Nature Biotechnology __34__, p 525--527 (2016).

Scripts reproducing all the results of the paper are available [here](https://github.com/pachterlab/kallisto_paper_analysis).

__kallisto__ quantified bulk RNA-Seq can be analyzed with [__sleuth__](https://github.com/pachterlab/sleuth/).

__kallisto__ can be used together with [__bustools__](https://bustools.github.io/) to pre-process single-cell RNA-seq data. See the [kallistobus.tools](https://www.kallistobus.tools/) website for instructions.

## Manual

Please visit http://pachterlab.github.io/kallisto/manual.html for the manual.

## License

__kallisto__ is distributed under the BSD-2 license. The license is distributed with __kallisto__ in the file `license.txt`, which is also viewable [here](https://pachterlab.github.io/kallisto/download). Please read the license before using __kallisto__.

## Getting help

For help running __kallisto__, please post to the [kallisto-and-applications Google Group](https://groups.google.com/forum/#!forum/kallisto-and-applications).

## Reporting bugs

Please report bugs to the [Github issues page](https://github.com/pachterlab/kallisto/issues)

## Development and pull requests

We typically develop on separate branches, then merge into devel once features
have been sufficiently tested. `devel` is the latest, stable, development
branch. `master` is used only for official releases and is considered to be
stable. If you submit a pull request (thanks!) please make sure to request to
merge into `devel` and NOT `master`. Merges usually only go into `master`, but
not out.
