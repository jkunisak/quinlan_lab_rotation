#!/bin/bash

## Assign variables
bed_file_dir="$1"
pathoscore_script="$2"
pathoscore_clinvar_reference="$3"
outputdir="$4"

## Print out the variables 
echo "bed file: " $bed_file_dir"/new_CCR.bed.gz"
echo "pathoscore script path: " $2
echo "pathoscore clinvar reference truth set directory: " $3

rm -rf $bed_file_dir/new_CCR.bed.gz*
rm -rf $outputdir/*

## Zip the bed file 
echo "gzipping the bed file"
bgzip $bed_file_dir/new_CCR.bed

## Index the bed file 
echo "index gzipped bed file"
tabix -p bed $bed_file_dir/new_CCR.bed.gz

## Annotate benign 
echo "annotate benign pathoscore"
python $pathoscore_script annotate --prefix benign --scores $bed_file_dir/new_CCR.bed.gz:new_CCR:6:max $pathoscore_clinvar_reference/clinvar-benign.20170905.vcf.gz

## Annotate pathogenic 
echo "annotate pathogenic pathoscore"
python $pathoscore_script annotate --prefix pathogenic --scores $bed_file_dir/new_CCR.bed.gz:new_CCR:6:max $pathoscore_clinvar_reference/clinvar-pathogenic.20170905.vcf.gz

## Run pathoscore evaluate 
echo "run pathoscore evaluate"
#python $pathoscore_script evaluate $pathoscore_output_dir/pathogenic.vcf.gz $pathoscore_output_dir/benign.vcf.gz
python $pathoscore_script evaluate -s new_CCR pathogenic.vcf.gz benign.vcf.gz

mv pathoscore.*.* $outputdir
mv pathoscore.csv $outputdir
mv pathoscore.html $outputdir
mv pathogenic.* $outputdir
mv benign.* $outputdir

#cat <(head -n 1 ccrs.autosomes.v2.20180420.bed) <(awk '$1 == 20' ~/Desktop/ccrs.autosomes.v2.20180420.bed) > ~/Desktop/quinlan_lab_rotation/chr20/chr20_jim_ccr.txt
