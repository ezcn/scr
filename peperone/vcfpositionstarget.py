SpcDct={}
for line in open("data/id_list.txt", "r" ) : 
	x=line.rstrip().split("\t")
	#print x  
	SpcDct[x[9]]=x[8]

#print SpcDct

count=0 
for line in open ("names.merged" , "r" ) : 
	x=line.rstrip()
	if x in SpcDct: 
		resline=[ count , x, SpcDct[x]]
	else: resline=[ count , x, "failed"]
	print "\t".join(map(str, resline) )
	count +=1

