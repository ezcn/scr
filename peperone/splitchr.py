import gzip 
import re 
import sys

rstart=int(sys.argv[1])
rend=int(sys.argv[2])

for n in range (rstart, rend): 
#	f=gzip.open("/home/enza/Desktop/peperone/data/all.mergedSNPs.vcf.gz", "r" )  
#	out=open("/home/enza/Desktop/peperone/data/%s.all.mergedSNPs.vcf.gz" %(n), "w")
	f=open("/home/enza/oogaprotocol/PEPERONE/geneticmap/genetic.map", "r" )  
	out=open("/home/enza/oogaprotocol/PEPERONE/geneticmap/%s.genetic.map" %(n), "w")
	sys.stdout=out 
	for line in f: 
#		if re.match("#", line): print line.rstrip() 
		if re.match("Chr", line): pass
		else: 
			x=line.split() 
			if int(x[0])==n: print line.rstrip()  
	f.close() 
