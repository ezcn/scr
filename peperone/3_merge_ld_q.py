import re 
d_info={}
for line in open ("/home/enza/oogaprotocol/PEPERONE/pop4/list.complete.pop4") :
	x=line.rstrip().split() 
	species=x[1]; clade=x[4]; mytype=x[5] 
	if not species =="annuum": 
		if not species in d_info: d_info[species]=[clade, mytype ]

d_info["annuum.sf"]=["annum","domesticated" ]
d_info["annuum.bf"]=["annum","domesticated" ]

with open("/home/enza/oogaprotocol/PEPERONE/pop4/species.q") as f:
    specieslist= f.read().splitlines()

title=["CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B", "R2", "species", "clade", "type" ]
print "\t".join(map(str, title ))
for sp in specieslist:
	if not (sp=="eximium" or sp=="galapagoense") :
		filename= "/home/enza/oogaprotocol/PEPERONE/6_ld_species_q/%s.ld" % (sp)
		for sline in open (filename): 
			y=sline.rstrip().split()
			if not re.search("CHR", sline): 
				res=y+[sp]+d_info[sp]
				print "\t".join(map(str, res )) 

		  
