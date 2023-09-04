# ornaments
__Ornaments__ is a wrapper program around kallisto to facilitate
variant-aware pseudoalignment of RNA-Seq data to obtain allele-specific read counts 
or perform allele specific transcript quantification. __ornaments__ takes as input
a list of variants (both SNPs and indels are supported) and a reference
transcriptome to create a modified reference that can accommodate allele-specific
alignments. __ornaments__ requires that the variants
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
Example transcriptome files and gene annotation files necessary for the scripts can be found at: [here](https://www.gencodegenes.org/human/).

You can download a list of variant files [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/). The files of interest are the .vcf.gz files which contain a list of all the phased variants, separated by chromosome.

For convenience, we have included a set of example files in the _examples_ folder, which we will use in the workflow examples below.

## Workflow
In order to run __ornaments__, you need to have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed. Once you have installed conda, install the [biopython](https://anaconda.org/conda-forge/biopython) and [gtfparse](https://anaconda.org/bioconda/gtfparse) packages:

```
conda create --name ornaments
conda activate ornaments
conda install -c conda-forge biopython
conda install -c bioconda gtfparse
```

To begin, we must translate our set of genomic variants into transcriptomic variants. 
The provided script `convert_genomic_vcf_to_transcriptomic_vcf.py` will transform
genomic variants (one chromosome at a time, in vcf format) to transcriptomic variants.

Example Usage:

```
python3 convert_genomic_vcf_to_transcriptomic_vcf.py \
-v examples/truncated_chr1.vcf.gz \
-g examples/truncated_v36.annotation.gtf \
-t examples/truncated_transcriptome.fasta \
-o examples/truncated_chr1.transcriptome.vcf \
-c examples/truncated_chr1.transcriptome.coords
```


Here, `examples/truncated_chr1.vcf.gz`is the example genome VCF file, `examples/truncated_v36.annotation.gtf` is the gene to
transcript annotation file, and `examples/truncated_transcriptome.fasta` is the transcriptome. 
The outputs of the script are `examples/truncated_chr1.transcriptome.vcf` which contains the remapped variants, 
and `examples/truncated_chr1.transcriptome.coords` which keeps track of coordinate conversions.

Please note that a compatible reference transcriptome and transcript annotation file should be used. 
If you are getting both files from the same source, make sure that the version is the same. If there 
are inconsistencies, this will be reflected in the output of the script. If there are no inconsistencies,
you should see the following output from the script:

```
Number of SNPs inconsistent with conversion:  0
```

Upon creating the transcriptomic variants, we need to generate a modified reference
transcriptome to use with the ornaments software. If your variants span more than one chromosome,
then we need to merge the transcriptome variant files from the previous step.
The `merge_transcriptome_vcfs.py` script is provided for this.

Example Usage:
```
python3 merge_transcriptome_vcfs.py -i in_files -o transcriptome.vcf -s sample_file
```

Here, `in_files` is a file containing the names of all the files to be merged, with one file name specified in each line.
In the examples folder, there is only one chromosome which we consider, so running this step is technically unnecessary.
`transcriptome.vcf` is the output file, and `sample_file` is a file with a list of sample names (separated by newlines)
which are considered, since only variants that are heterozygous in at least one of the samples will be retained. 
Alternatively, one can specify a single sample in place of a sample file.

Finally, we must also make sure that the output transcriptome variant file is sorted with respect to each transcript,
i.e. the list of variants for each transcript appears in sorted order. If the VCF file you started with is already sorted, 
then great, there is nothing more to do here! Otherwise, you can use this command to sort the transcriptome VCF, taken from 
[here](https://www.biostars.org/p/299659/):

```
cat transcriptome.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > sorted.transcriptome.vcf
```

You should replace `transcriptome.vcf` with your file names. Upon merging and sorting, we can use the provided script 
`create_personalized_transcriptome.py` to modify the reference.

Example Usage to Create Diploid Transcriptome:
```
python3 create_personalized_transcriptome.py \
-f examples/truncated_transcriptome.fasta \
-v examples/truncated_chr1.transcriptome.vcf \
-t diploid \
-o SAMPLE.diploid \
-s SAMPLE
```

Here, `diploid` denotes that a diploid transcriptome will be created. This will result in a file called `SAMPLE.diploid.fa`,
which contains left (suffixed by `_L`) and right (suffixed by `_R`) alleles for each transcript, wherein the phased variants
have been imputed into each transcript. The name `SAMPLE` that should be used is the name that refers to the sample in the VCF file.

Example Usage to Create Ornament Transcriptome:
```
python3 create_personalized_transcriptome.py \
-f examples/truncated_transcriptome.fasta \
-v examples/truncated_chr1.transcriptome.vcf \
-t ornament \
-o SAMPLE.ornament \
-s SAMPLE
```

Here, `diploid` denotes that a diploid transcriptome will be created. This will result in a file called `SAMPLE.diploid.fa`,
which contains left (suffixed by `_L`) and right (suffixed by `_R`) alleles for each transcript, wherein the phased variants
have been imputed into each transcript. The name `SAMPLE` that should be used is the name that refers to the sample in the VCF file.

Here, `ornaments` denotes that an ornament transcriptome will be created. This will result in a file called `SAMPLE.ornament.fa`, 
which can be used by __ornaments__. __Ornaments__ expands the functionality of __kallisto__, meaning it functions differently 
only if `SAMPLE.ornament.fa` is provided as input. Otherwise, with a reference transcriptome, it should behave the exact same as __kallisto__.

Running Ornaments:

```
build/src/ornaments index -i SAMPLE.ornament.fa SAMPLE.ornament.index
```

This will construct an ornament index file, called `SAMPLE.ornament.index`, which will be used for pseudoalignment in the quantification. 

To perform variant-aware quantification, we need to provide a few more parameters. To run with paired-end data that is specific to a given `SAMPLE`:

```
build/src/ornaments quant \
-i SAMPLE.ornament.index \
-o output_dir \
--vcf sorted.transcriptome.vcf \
--sample SAMPLE \
SAMPLE_1.fastq SAMPLE_2.fastq
```

Here, the `vcf` parameter is used to read the variant information, and the `sample` parameter tells 
__ornaments__ which sample to focus on for allele specificity, as different samples will have
varying patterns of heterozygosity. This will result in an output file `output_dir/allele_counts.txt`,
which contains the allele specific counts for each heterozygous variant in the input VCF file, and `output_dir/tpms.txt`,
which contain TPM estimates for each transcript in the input transcriptome.

To do the same with single-end RNA-seq data, we additionally need to specify the fragment length mean and standard deviation used during sequencing:

```
build/src/ quant \
-i SAMPLE.ornament.index \
-o output_dir \
--vcf sorted.transcriptome.vcf \
--sample SAMPLE \
--single \
-l 160 -s 25 \
SAMPLE.fastq
```

In this example the fragment length mean is 160 bp and the standard deviation is 25 bp (you can use these values for Geuvadis data).

Ornaments is built as an addition to __kallisto__. All credit for the __kallisto__ method is given to its original authors.