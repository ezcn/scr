import sys, re

infile = sys.argv[1] 

title=['id', 'clade', 'nat', 'con', 'species', 'type', 'miss', 'id2', "cluster", "cval"]
print  "\t".join(title )

data_admix={} 

countid=0
for line in open (infile, "r" ): 
	x=line.rstrip().split() 

	if not re.match("id", line) : 
		countcluster=1 ; countid+=1 
		for k in x[8:]: 
			resline = ["%s_%s" % (i, countid) for i in  x[0:8]]  + [ "c%s" %(countcluster) , k ]	
			countcluster+=1 
			print "\t".join( map (str, resline ) ) 

