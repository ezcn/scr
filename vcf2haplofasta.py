#!/usr/bin/python
import os
import sys

"""
extract haplotypes in fasta format  from vcf files
only SNPs retained
mylistofindividulas = single column list of individuals (should correspond to id in the vcf file ) 


usage : ./extracthaplotypes_bender.py 2:222-222 mycvf mylistofindividulas 

"""


myregion=sys.argv[1]
myvcf=sys.argv[2]
indlist=sys.argv[3]



cmdl='for i in $(cat %s); do for n in 1 2; do echo  ">"$i-$n"\n"$(tabix -h  %s  %s | vcfkeepsamples - $i | vcffilter -f "VT = SNP" | vcfgenotypes - | tr \":\" \"\\t\"  | cut -f7 | tr \"/\" \"\\t\"  | cut -f$n | perl -pi -e \"s/\\n//g\"); done; done ' % (indlist, myvcf , myregion  )

os.system(cmdl)
#print cmdl
