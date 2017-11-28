import re 



#types=["coltivato", "selvatico", "domestico"]
types= ["wild" , "domesticated"] 
for t in types: 
	for c in range (1, 13): 
		pathtodir="/home/enza/oogaprotocol/PEPERONE/3_ld"
#		pathtodir="/home/enza/oogaprotocol/PEPERONE/ldblocks"
		myfile="%s.%s.all.merged.ld" %(c, t )
#		myfile="%s.%s.all.merged.blocks.det" %(c, t )
		with open ("%s/%s" % (pathtodir, myfile) , "r" ) as f :
			next(f)
			for line in f:  
				x=line.split() 
				x.append(t)	
				print "\t".join(x) 
				 
