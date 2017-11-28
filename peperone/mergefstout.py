import re 



types=["selvatico.coltivato", "selvatico.domestico", "domestico.coltivato"]
for t in types: 
	for c in range (1, 13): 
		#pathtodir="/home/enza/oogaprotocol/PEPERONE/segmentfst"
		pathtodir="/home/enza/oogaprotocol/PEPERONE/fst_pruned"
#		myfile="%s.%s.fst.bed" %(c, t )
		myfile="%s.%s.pruned.fst" %(c, t )
		with open ("%s/%s" % (pathtodir, myfile) , "r" ) as f :
			#next(f)
			for line in f:  
				x=line.split() 
				x.append(t)	
				print "\t".join(x) 
				 
