# Ornaments
Ornaments is a lightweight modification of kallisto to facilitate
variant-aware pseudoalignment of RNA-Seq data to obtain allele-specific read counts 
or perform allele specific transcript quantification. Ornaments takes as input
a list of variants (both SNPs and indels are supported) and a reference
transcriptome to create a modified reference that can accommodate allele-specific
alignments. Ornaments requires that the variants
are listed in terms of transcriptomic coordinates. We provide helper scripts
to prepare the ornaments reference as well as transform the variant coordinates.

## Compiling
To compile ornaments, cmake is required. First, clone the repository and go into the
ornaments folder. 

```
git clone https://gitlab.com/aadduri/ornaments.git
cd ornaments
```

Next, we need to create a build folder, which will contain the __ornaments__ executable.

```
mkdir build
cd build
```

Finally, we will build the executable using cmake. If you do not have cmake installed, you can get it by running 

```
brew install cmake
```

Then, run:

```
cmake ..
make
```

An executable should be built with the name _ornaments/build/src/ornaments_. If you would rather install the software
in your home directory, rather than the above, you should run:
```
cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$HOME
make install
```

which will build an executable at _$HOME/bin/ornaments_.

## Data
Example transcriptome files and gene annotation files necessary for the scripts can be found [here](https://www.gencodegenes.org/human/).

You can download a list of variant files [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/). The files of interest are the .vcf.gz files which contain a list of all the phased variants, separated by chromosome.

For convenience, we have included a set of example files in the _examples_ folder, which we will use in the workflow examples below.

## Setting up the dependencies
In order to run __Ornaments__, you need to have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed. Once you have installed conda, install the [biopython](https://anaconda.org/conda-forge/biopython) and [gtfparse](https://anaconda.org/bioconda/gtfparse) packages:

```
conda create --name ornaments
conda activate ornaments
conda install -c conda-forge biopython
conda install -c bioconda gtfparse
```

## Preparation of variants

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
The script `merge_transcriptome_vcfs.py` merges the transcriptome variant files from the previous step,
if your variants span more than one chromosome.

```
Usage: python3 merge_transcriptome_vcfs.py -i <file_list> -o <file_vcf> -s <file_sample>

-i <file_list>: 	input file listing the names of the files to be merged, with one file name in each line
-o <file_vcf>: 		output file `transcriptome.vcf` 
-s <file_sample>: 	input file with a list of sample names, one sample name in each line
```

Only variants that are heterozygous in at least one of the samples will be retained. We also provide an alternative
usage for this script, where you can specify a single sample instead of a file containing sample names:

```
Usage: python3 merge_transcriptome_vcfs.py -i <file_list> -o <file_vcf> -s <sample_name>

-i <file_list>: 	input file listing the names of the files to be merged, with one file name in each line
-o <file_vcf>: 		output file `transcriptome.vcf` 
-s <sample_name>: 	the sample name in the VCF file, i.e. "HG00405"
```

In the examples folder, there is only one chromosome which we consider, so running this step is not necessary.

### Sorting (if variants are not sorted in transcriptome variant files)
If the VCF file you started with is already sorted, you can skip this part. 
Otherwise, you can use this command to sort the variants for each transcript 
in the transcriptome variant file, to sort the transcriptome VCF, taken from 
[here](https://www.biostars.org/p/299659/).

```
cat transcriptome.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > sorted.transcriptome.vcf
```

You should replace `transcriptome.vcf` with your file names. Upon merging and sorting, we can use the provided script 

## Ornament personalized transcriptome
The script `create_personalized_transcriptome.py` modifies the reference transcriptome and produces ornament personalized transcriptome or diploid personalized transcriptome.

```
Usage: python3 create_personalized_transcriptome.py -f <file_tr> -v <file_vcf> -t <opt> -o <file_ptr> -s <SAMPLE>

-f <file_tr>: 		input reference transcriptome file
-v <file_vcf>:		input VCF file for variants in transcriptome coordinates
-t <opt>: 		set opt to ornaments for ornament personalized transcriptome or 
						to diploid for diploid personalized transcriptome
-o <file_ptr>:		output file containing ornament or diploid personalized transcriptome
-s <SAMPLE>: 		the sample name in the VCF file, i.e. "HG00405"
-k <kmer_len>: 		the length of kmer to use for ornament transcriptome construction, default = 31
```

Diploid personalized transcriptome contains left (suffixed by `_L`) and right 
(suffixed by `_R`) alleles for each transcript, 
wherein the phased variants have been imputed into each transcript. 

`SAMPLE.ornament.fa`, output of this script, input file to Ornaments 
__Ornaments__ expands the functionality of kallisto, meaning it functions differently 
only if `SAMPLE.ornament.fa` is provided as input. Otherwise, with a reference transcriptome, 
it should behave the exact same as kallisto.

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

### index


```
Usage: build/src/ornaments index -i <out_index> -k <k> <in_fasta>

Arguments:
-i <out_index>:		output file that contains an ornament index 
<in_fasta>:		    input file for personalized transcriptome
-k <k>: 		    k-mer length to be used for index construction, default = 31

```
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
--sample <SAMPLE>: 	sample name in the population data for quantification, i.e. "HG00405"
<file_fastq1>:		FASTQ file 1 in paired-end reads 
<file_fastq2>: 		FASTQ file 2 in paired-end reads
--single:		flag for single-end reads 
-l <d_l>: 		fragment length mean for single-end reads
-s <d_s>: 		standard deviation for single-end reads
<file_fastq>:		FASTQ file for single-end reads
```

Example usage to quantify with paired-end reads for a given `SAMPLE`:

```
build/src/ornaments quant -i examples/HG00405.ornament.index -o output_dir \
--vcf examples/sorted.transcriptome.vcf --sample HG00405 examples/SAMPLE_1.fastq examples/SAMPLE_2.fastq
```

Example usage to quantify with single-end reads for a given `SAMPLE`: 
```
build/src/ornaments quant -i examples/HG00405.ornament.index -o output_dir \
--vcf examples/sorted.transcriptome.vcf --sample HG00405 --single -l 150 -s 25 examples/SAMPLE.fastq
```

[comment]: When I run build/src/ornaments without any options, I get kallisto's usage. Once we finalize this document, this should be replaced with the Usage for Ornaments.

Ornaments is built as an addition to kallisto. All credit for the kallisto software is given to its original authors.
