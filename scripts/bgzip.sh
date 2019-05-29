#!/bin/bash

## Assign variables 
bed_file_dir="$1" 

## Clean out the directory  
rm -rf $bed_file_dir/*.gz*
rm -rf $bed_file_dir/storage/*

## First copy the bed files to a storage directory 
cp $bed_file_dir/*.bed $bed_file_dir/storage

## bgzip the bed files 
ls $bed_file_dir/*.bed | awk '{print "bgzip "$1""}' | bash 

## Index the zipped bed files 
ls $bed_file_dir/*.bed.gz | awk '{print "tabix -p bed "$1""}' | bash

## Move the original bed files back to the working directory 
cp $bed_file_dir/storage/* $bed_file_dir 
