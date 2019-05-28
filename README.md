# quinlan_lab_rotation
This is the repository that contains the data and files needed to generate the 2nd version of CCRs. This version includes variants from gnomad exome and whole genome sequencing samples. 

The goal of this project is to generate a new model for CCRs that will account for different SIFT/PolyPhen scores and VAF. This will allow for the new CCR model to identify regions enriched in Clinvar pathogenic variants. 

Required files include: 
1) exter.py function to read the gff file 
2) pathoscore.py function from [PathoScore repository](https://github.com/quinlan-lab/pathoscore/blob/master/pathoscore.py)
3) PathoScore Clinvar truth sets [see here to make the Clinvary truth-set directory](https://github.com/quinlan-lab/pathoscore/tree/master/truth-sets/GRCh37/clinvar)
4) 
