# quinlan_lab_rotation
This is the repository that contains the data and files needed to generate the 2nd version of CCRs. This version includes variants from gnomad exome and whole genome sequencing samples. 

The goal of this project is to generate a new model for CCRs that will account for different SIFT/PolyPhen scores and VAF. This will allow for the new CCR model to identify regions enriched in Clinvar pathogenic variants. 

**Required files include:**
1) new_CCR_jason.ipynb = jupyter notebook for generation of new CCR models 
2) exter.py function to read the gff file 
3) pathoscore.py function from [PathoScore repository](https://github.com/quinlan-lab/pathoscore/blob/master/pathoscore.py)
4) PathoScore Clinvar truth sets: [see here to make the Clinvary truth-set directory](https://github.com/quinlan-lab/pathoscore/tree/master/truth-sets/GRCh37/clinvar)

**Output files:**
1) pathoscore results will be provided in the pathoscore_results output directory 
2) new_CCR.bed.gz = bed file with the CCR windows 
3) new_CCR.bed.gz.tbi = indexed bed file 
4) new_CCR.txt = txt file with the CCR windows 
