# Ornaments
Ornaments is a lightweight modification of kallisto that facilitates
variant-aware pseudoalignment of RNA-Seq reads to obtain expected allele-specific read counts 
at heterozygous variant loci in addition to transcriptome quantification. 

Ornaments takes as input an ornament personalized transcriptome and reads.
Ornaments requires that the variants are listed in terms of transcriptomic coordinates. 
We provide helper scripts to transform the variant coordinates from genomic to transcriptomic coordinates
and to generate an ornament personalized transcriptome from 
a list of variants (both SNPs and indels are supported) and a reference transcriptome.

## Compiling
To compile ornaments, cmake is required. First, clone the repository and go into the
ornaments folder. 

```
git clone https://gitlab.com/aadduri/ornaments.git
cd ornaments
```

Next, we need to create a build folder. This folder will contain the ornaments executable.

```
mkdir build
cd build
```

Finally, we will build the executable using cmake. If you do not have cmake installed, you can get it by running 

```
brew install cmake
```

Then, run

```
cmake ..
make
```

An executable will be built with the path/name _ornaments/build/src/ornaments_. If you would rather install the software
in your home directory, rather than the above, you should run
```
cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$HOME
make install
```

This will build an executable at _$HOME/bin/ornaments_.

## Data
Example transcriptome files and gene annotation files necessary for the helper scripts can be found [here](https://www.gencodegenes.org/human/).

You can download a list of variant files [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/). The files of interest are the .vcf.gz files that contain a list of all the phased variants, separated by chromosome.

For convenience, we have included a set of example files in the _examples_ folder, which we will use in the workflow examples below.

## Setting up the dependencies
In order to run Ornaments, you need to have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed. Once you have installed conda, install the [biopython](https://anaconda.org/conda-forge/biopython) and [gtfparse](https://anaconda.org/bioconda/gtfparse) packages:

```
conda create --name ornaments
conda activate ornaments
conda install -c conda-forge biopython
conda install -c bioconda gtfparse
```

These packages are not needed to run _ornaments_ itself, but are needed for the helper scripts described below to prepare the inputs to _ornaments_.

## Preparation of variants


Since Ornaments requires that the variants are listed in terms of transcriptomic coordinates,
are ordered within each transcript, and are merged into a single file if stored across multiple files, 
we provide helper scripts to transform the variant coordinates from genomic to transcriptomic coordinates,
to sort the variants, and to merge the variant files. The output variant file from these scripts will be
used as input to construct an ornament personalized transcriptome.



### Converting genomic coordinates to transcriptome coordinates for variants

The script `convert_genomic_vcf_to_transcriptomic_vcf.py` transforms
the genome coordinates of variants to transcriptomic coordinates.


```
Usage: python3 convert_genomic_vcf_to_transcriptomic_vcf.py -v <file_vcf> -g <file_annot> -t <file_tr> -o <file_tr_vcf> -c <file_conv>

-v <file_vcf>:		input gzipped VCF file for variants in genome coordinates
-g <file_annot>: 	input gene-to-transcript annotation file
-t <file_tr>: 		input reference transcriptome file
-o <file_tr_vcf>: 	output VCF file for variants in transcriptome coordinates 
-c <file_conv>: 	output file with information on coordinate conversions
```

Example usage:

```
python3 convert_genomic_vcf_to_transcriptomic_vcf.py \
-v examples/truncated_chr1.vcf.gz \
-g examples/truncated_v36.annotation.gtf \
-t examples/truncated_transcriptome.fasta \
-o examples/truncated_chr1.transcriptome.vcf \
-c examples/truncated_chr1.transcriptome.coords
```


A compatible reference transcriptome and transcript annotation file should be used. 
If you are getting both files from the same source, make sure that the version is the same. If the
versions are the same, then you should see the following output from the script:


```
Number of SNPs inconsistent with conversion:  0
```

### Merging (if there are multiple transcriptome variant files)
If your variants are stored in multiple files (e.g., your variants span more than one chromosome and 
each file stores one chromosome at a time), these files should be merged into a single file.
The script `merge_transcriptome_vcfs.py` merges the transcriptome variant files from the previous step.


```
Usage: python3 merge_transcriptome_vcfs.py -i <file_list> -o <file_vcf> -s <file_sample>

-i <file_list>: 	input file listing the names of the files to be merged, with one file name in each line
-o <file_vcf>: 		output merged VCF file for variants in transcriptome coordinates
-s <file_sample>: 	input file with a list of sample names, one sample name in each line, if processing multiple samples,
                            or the sample name in the VCF file, if processing one sample
```

This script expects all files in `<file_list>` to be in VCF format, and each file should have data for all samples in `<file_sample>`.

Only variants that are heterozygous in at least one of the samples in `<file_sample>` will be retained. In the examples folder, we consider only one chromosome, but we should still run this script to remove variants which are not heterozygous in at least one sample.

Example usage:
```
echo "examples/truncated_chr1.transcriptome.vcf" > file_list.txt
python merge_transcriptome_vcfs.py -i file_list.txt -o examples/transcriptome.vcf -s HG00405
```

### Sorting (if variants are not sorted in transcriptome variant files)
If variants in the VCF file you started with are already sorted, you can skip this part. 
Otherwise, you can use this command to sort the variants for each transcript 
in the transcriptome variant file, to sort the transcriptome VCF, taken from 
[here](https://www.biostars.org/p/299659/).

```
cat examples/transcriptome.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > examples/sorted.transcriptome.vcf
```

You should replace `transcriptome.vcf` with your file name. 

## Ornament personalized transcriptome
The script `create_personalized_transcriptome.py` modifies the reference transcriptome with the variant information that was prepared above, and produces ornament personalized transcriptome. This script can also be used to generate a diploid personalized transcriptome.

```
Usage: python3 create_personalized_transcriptome.py -f <file_tr> -v <file_vcf> -t <opt> -o <file_ptr> -s <SAMPLE>

-f <file_tr>: 		input reference transcriptome file
-v <file_vcf>:		input VCF file for variants in transcriptome coordinates
-t <opt>: 		set <opt> to ornaments for ornament personalized transcriptome or 
                            to diploid for diploid personalized transcriptome
-o <file_ptr>:		output file containing ornament or diploid personalized transcriptome
-s <SAMPLE>: 		one desired sample from the VCF file, i.e. "HG00405" which should be contained in <file_vcf>
-k <kmer_len>: 		the length of kmer to use for ornament transcriptome construction, default = 31
```

Diploid personalized transcriptome contains left (suffixed by `_L`) and right 
(suffixed by `_R`) alleles for each transcript, 
wherein the phased variants have been imputed into each transcript. 


The output file of this script `<file_ptr>` contains the ornament personalized transcriptome to be used as an input file to Ornaments. 
Ornaments expands the functionality of kallisto, meaning it functions differently 
only if `<file_ptr>` is provided as input. Otherwise, with a reference transcriptome, 
it should behave exactly the same as kallisto.

Example usage to create ornament personalized transcriptome:
```
python3 create_personalized_transcriptome.py \
-f examples/truncated_transcriptome.fasta \
-v examples/truncated_chr1.transcriptome.vcf \
-t ornament \
-o SAMPLE.ornament.fa \
-s SAMPLE
```


Example usage to create diploid personalized transcriptome:
```
python3 create_personalized_transcriptome.py \
-f examples/truncated_transcriptome.fasta \
-v examples/truncated_chr1.transcriptome.vcf \
-t diploid \
-o SAMPLE.diploid.fa \
-s SAMPLE
```


## Ornaments 

Given the ornaments personalized transcriptome above, Ornaments is run with two successive commands: `index` for generating an ornaments index and `quant` for quantification of transcriptome and allele-specific expression levels at heterozygous variant loci.

### index


```
Usage: build/src/ornaments index -i <file_index> -k <k> <file_fasta>

Arguments:
-i <file_index>:	output file that contains an ornament index 
<file_fasta>:		input file for personalized transcriptome (ornament personalized transcriptome to run ornaments, normal reference transcriptome to run kallisto)
-k <k>:			k-mer length to be used for index construction, default = 31

```

Note that ornaments index can automatically detect whether or not the input `<file_fasta>` is an ornament personalized transcriptome or a normal reference transcriptome. If the input if a normal reference transcriptome, then kallisto will be run instead.

Example usage to construct an ornament index:

```
build/src/ornaments index -i SAMPLE.ornament.index SAMPLE.ornament.fa
```

### quant

```
Usage:  
build/src/ornaments quant -i <file_orn_ind> -o <dir_out> --vcf <file_vcf> --sample <SAMPLE> <file_fastq1> <file_fastq2> 
build/src/ornaments quant -i <file_orn_ind> -o <dir_out> --vcf <file_vcf> --sample <SAMPLE> --single -l <d_l> -s <d_s> <file_fastq>

Arguments:
-i <file_orn_ind>:	input file for ornament index 
-o <dir_out>: 		two output files will be added to the <dir_out> folder 
                        `allele_counts.txt` for expected allele specific-read counts 
                        `tpms.txt` for TPM estimates for transcripts 
--vcf <file_vcf>:	sorted.transcriptome.vcf the variant information, 
--sample <SAMPLE>: 	one desired sample from the VCF file, i.e. "HG00405" which should be contained in <file_vcf>
<file_fastq1>:		FASTQ file 1 in paired-end reads 
<file_fastq2>: 		FASTQ file 2 in paired-end reads
--single:		flag for single-end reads 
-l <d_l>: 		fragment length mean for single-end reads
-s <d_s>: 		standard deviation for single-end reads
<file_fastq>:		FASTQ file for single-end reads
```

Example usage to quantify with paired-end reads for a given `SAMPLE`:

```
build/src/ornaments quant \
-i examples/HG00405.ornament.index \
-o output_dir \
--vcf examples/sorted.transcriptome.vcf \
--sample HG00405 \
examples/SAMPLE_1.fastq examples/SAMPLE_2.fastq
```

Example usage to quantify with single-end reads for a given `SAMPLE`: 
```
build/src/ornaments quant \
-i examples/HG00405.ornament.index \
-o output_dir \
--vcf examples/sorted.transcriptome.vcf \
--sample HG00405 \
--single \
-l 150 \
-s 25 \
examples/SAMPLE.fastq
```

Ornaments is built as an addition to kallisto. All credit for the kallisto software is given to its original authors.
