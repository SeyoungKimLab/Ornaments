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

Finally, we will build the executable using cmake.
```
cmake ..
make
```

An executable should be built with the name _ornaments/build/src/ornaments_.

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
the genome coordinates of variants to transcriptomic coordinates
(one chromosome at a time, in vcf format). 

[comment]: Above, is it necessary to say one chromosome at a time? can we delete this?


```
Usage: python3 convert_genomic_vcf_to_transcriptomic_vcf.py -v <file_vcf> -g <file_annot> -t <file_tr> -o <file_tr_vcf> -c <file_conv>

-v <file_vcf>:		input VCF file for variants in genome coordinates
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


[Comment]: In the example above, if the VCF file has to be gzipped, then say that in the Usage.

A compatible reference transcriptome and transcript annotation file should be used. 
If you are getting both files from the same source, make sure that the version is the same. If there 
are inconsistencies, this will be reflected in the output of the script. If there are no inconsistencies,
you should see the following output from the script:

[Comment]: the sentence above [If there are inconsistencies, this will be reflected ....] Could you revise this tfor clarity? i.e., what is the output of the script here, remapped or coord file above? what do you mean by [reflected]?

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

[Comment]: in the code block above, the description of file_sample is unclear. What are sample names? Are they filenames, or sample names inside a file? And, in the sentence below, [one can specify a single sample in place of a sample file], this is not clear. Is this a filename for that single sample? What is the format of that file?


Only variants that are heterozygous in at least one of the samples will be retained. 
Alternatively, one can specify a single sample in place of a sample file.



Example usage:
```
python3 merge_transcriptome_vcfs.py -i in_files -o transcriptome.vcf -s sample_file
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
-s <SAMPLE>: 		the name that refers to the sample in the VCF file.
```

[Comment]: In the code block above, the description of SAMPLE is not clear. See my comments on SAMPLE above as well. They are related comments.

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
-o SAMPLE.ornament \
-s SAMPLE
```


Example usage to create diploid personalized transcriptome:
```
python3 create_personalized_transcriptome.py \
-f examples/truncated_transcriptome.fasta \
-v examples/truncated_chr1.transcriptome.vcf \
-t diploid \
-o SAMPLE.diploid \
-s SAMPLE
```

[Comment]: SAMPLE.ornament and SAMPLE.diploid in the example above. Shouldn't this be SAMPLE.ornament.fa etc?


## Ornaments 

### index

[Comment]: Does `ornaments index' detect whether it is getting ornament or diploid personalized transcriptome automatically? Or does it have to be specified in the file name like *.ornament and *.diploid?

```
Usage: build/src/ornaments index -i <file_orn> <file_orn_ind>

Arguments:
-i <file_orn>:		input file for personalized transcriptome
<file_orn_ind>:		output file that contains an ornament index 

```
Example usage to construct an ornament index:

```
build/src/ornaments index -i SAMPLE.ornament.fa SAMPLE.ornament.index
```

[Comment]: I see in kallisto's install.md 1. Execute cmake. There are a few options: - `-DCMAKE_INSTALL_PREFIX:PATH=$HOME` which will put kallisto in `$HOME/bin` as opposed to the default (`/usr/local/bin`). Do the same and use ornaments in the usage and example instead of build/src/ornaments. Make sure you try it and this works.

[Comment]: I don't see an example fastq file supplied in the example folder. Does kallisto have an example?


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
--sample <SAMPLE>: 	sample ID in the population data for quantification 
<file_fastq1>:		FASTQ file 1 in paired-end reads 
<file_fastq2>: 		FASTQ file 2 in paired-end reads
--single:		flag for single-end reads 
-l <d_l>: 		fragment length mean for single-end reads
-s <d_s>: 		standard deviation for single-end reads
<file_fastq>:		FASTQ file for single-end reads
```

[Comment]: Again, the role and format of SAMPLE is unclear.

Example usage to quantify with paired-end reads for a given `SAMPLE`:

```
build/src/ornaments quant -i ornament_index -o output_dir \
--vcf sorted.transcriptome.vcf --sample SAMPLE SAMPLE_1.fastq SAMPLE_2.fastq
```

Example usage to quantify with single-end reads for a given `SAMPLE`: 
```
build/src/ornaments quant -i ornament_index -o output_dir \
--vcf sorted.transcriptome.vcf --sample SAMPLE --single -l 160 -s 25 SAMPLE.fastq
```

[comment]: About [the fragment length mean 160 bp and the standard deviation 25 bp can be used for the Geuvadis data], is the example file a Geuvadis sample? I don't think it is? I think the mean length and standard deviation come from each dataset, in which case this is dataset-specific and the numbers from that dataset should be used. If this is the case, you can leave out this sentence.

[comment]: I don't see SAMPLE_1.fastq, SAMPLE_2.fastq, and SAMPLE.fastq in the examples folder. The folder should include these files. One should be able to copy and paste the line above to the command line and hit enter.

[comment]: When I run build/src/ornaments without any options, I get kallisto's usage. Once we finalize this document, this should be replaced with the Usage for Ornaments.

[comment]: kallisto has lots of other options. Does Ornaments inherit this? Some are pretty basic, like k in k-mer and definitely should be included.


Ornaments is built as an addition to kallisto. All credit for the kallisto software is given to its original authors.
