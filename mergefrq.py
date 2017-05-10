import sys
import re

files=sys.argv[1:]

# merges the output of vcftools --freq (or --counts) 
#only for biallelic loci 
#as many pop as you want 
#out: merged.freq.out (merged allele frequencies) 
#out: excluded.freq.out  (not same refernce allele or site not in all files)  

#read files and store loci in site_dz
site_dz={}
for fname in files:
	with open(fname) as f:
   		next(f)
		for line in f: 	
			x=line.rstrip().split()
			if (x[0], x[1]) in site_dz:
               			site_dz[(x[0], x[1])].append(x[4:]) 
      	          	else:
                        	site_dz[(x[0], x[1])]= [x[4:]]

#
# 
excludedloci=[] # loci that are not present in all files

out=open("merged.freq.out", "w")  
sys.stdout=out   

###
title=["chr", "loc", "ref", "alt"]
countpop=1
for n in range(len(files)): 
	title+=["ref_p%s" %(countpop) , "alt_p%s" %(countpop) ]
	countpop+=1

print  "\t".join(title)

###
for loc in site_dz:
	#print loc ,site_dz[loc],  len(site_dz[loc]), len(files) 
	if len(site_dz[loc])<len(files):
		reason=loc+ ("not in all files ",)
		excludedloci.append(reason)

	else:   
		#print loc,  site_dz[loc]
		res=[l for l in loc]
		for  info in  site_dz[loc]: 
			if not len(info)==2: 
				info.append("X:0") 
		#print site_dz[loc] 

		#check alleles are the same 
		allelearethesame=True 
		########## TO FIX ###################
		#for n in len(files): 
			#zipped=zip(site_dz[loc][0], site_dz[loc][1], site_dz[loc][2]) 
		zipped=zip(site_dz[loc][0], site_dz[loc][1]) 
		for allele in zipped:
			#print "****", allele
			theallele=[aa.split(":")[0] for aa in allele if not aa.split(":")[0]=="X"]	 
			#print theallele 
			if not len(set(theallele)) <= 1: allelearethesame=False 
		if  allelearethesame: 
			refall= zipped[0][0].split(":")[0]
                        altall= zipped[1][0].split(":")[0]
                        res+=[refall, altall]
			ref=[aa.split(":")[1] for aa in zipped[0]] #print "yay" , zipped 
			alt=[aa.split(":")[1] for aa in zipped[1]]
			refalt_bypop= [item for sublist in zip(ref,alt) for item in sublist]
			res+=refalt_bypop
			print "\t".join(res)
		else: 
			reason=loc+ ("not same allele", ) 
			excludedloci.append(reason)

out.close() 

out=open("excluded.freq.out", "w" ) 
sys.stdout=out 
for ex in excludedloci: 
	res=[i for i in ex]
	print "\t".join(res)  

