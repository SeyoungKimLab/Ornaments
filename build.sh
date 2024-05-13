#!/bin/bash

mkdir build
cd build
cmake ..
make

# Copy binaries and scripts to the Conda package bin directory
cp src/ornaments "${PREFIX}/bin/"
cp ../merge_transcriptome_vcfs.py "${PREFIX}/bin/merge_transcriptome_vcfs"
cp ../convert_genomic_vcf_to_transcriptomic_vcf.py "${PREFIX}/bin/convert_genomic_vcf_to_transcriptomic_vcf"
cp ../create_personalized_transcriptome.py "${PREFIX}/bin/create_personalized_transcriptome"

# Make the python scripts executable
chmod +x "${PREFIX}/bin/merge_transcriptome_vcfs"
chmod +x "${PREFIX}/bin/convert_genomic_vcf_to_transcriptomic_vcf"
chmod +x "${PREFIX}/bin/create_personalized_transcriptome"
