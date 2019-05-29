#!/bin/bash

## Assign variables 
outputdir="$1"

mv pathoscore.*.* $outputdir
mv pathoscore.csv $outputdir
mv pathoscore.html $outputdir
mv pathogenic.* $outputdir
mv benign.* $outputdir
