
d_info={}
for line in open ("/home/enza/oogaprotocol/PEPERONE/pop4/list.complete.pop4") :
	x=line.rstrip().split() 
	myid=x[0]; species=x[1]; clade=x[4]; mytype=x[5] 
 	d_info[myid]=[species,clade, mytype ]

with open("/home/enza/oogaprotocol/PEPERONE/pop4/species") as f:
    specieslist= f.read().splitlines()

title=["FID1", "IID1",  "FID2",  "IID2",  "RT", "EZ",  "Z0", "Z1", "Z2", "PI_HAT", "PHE", "DST", "PPC", "RATIO", "species", "clade", "type"]
#title=["FID",  "IID", "OHOM", "EHOM", "NHOM", "F", "species", "clade", "type" ]
print "\t".join(map(str, title ))
for sp in specieslist:
	if not (sp=="eximium" or sp=="galapagoense") :  
		#filename= "/home/enza/oogaprotocol/PEPERONE/8_inbreeding/het/%s.het" % (sp)
		#filename= "/home/enza/oogaprotocol/PEPERONE/8_inbreeding/ibs/%s.genome" % (sp)
		filename= "/home/enza/oogaprotocol/PEPERONE/8_inbreeding/ibs_nonfounder/%s.genome" % (sp)
		for sline in open (filename): 
			y=sline.rstrip().split()
			sid=y[0]
			if sid in d_info: 
				res=y+d_info[sid]
				print "\t".join(map(str, res )) 

		  
