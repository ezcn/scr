import os.path

import re 


speciestoexclude=["annuumvarglabrisculum", "galapagoense", "eximium", "praetermissum", "baccatumvarbaccatum"]

title=["CHROM","POS", "WEIR_AND_COCKERHAM_FST", "pair"]
print "\t".join(map(str, title))

for line in open ("/home/enza/oogaprotocol/PEPERONE/pop4/species.pairs") : 
	x=line.rstrip().split() 
	if not x[0] in speciestoexclude: 
		if not x[1] in speciestoexclude:  
			filename="/home/enza/oogaprotocol/PEPERONE/6_fst_species/%s.%s.weir.fst" % (x[0], x[1])
			if os.path.isfile(filename): 
				with open(filename) as f:
		    			next(f)
		    			for yline in f: 
						if not re.search("-nan", yline): 
							y=yline.rstrip().split() 
							res=y+["%s-%s" % (x[0], x[1]) ] 
							print "\t".join(map(str, res ))
