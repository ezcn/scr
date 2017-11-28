import gzip, sys , re 

infile=sys.argv[1] 

map_data={}
for line in gzip.open ("/home/enza/oogaprotocol/PEPERONE/5_neld/1kg/chr18.interpolated_genetic_map.gz", "r" ): 
	y=line.rstrip().split() 
	map_data[y[0]]=y[2]


for rline in open(infile, 'r') : 
	x=rline.rstrip().split()
	if re.match("CHR", rline) : 
		x+=['pos_a', 'pos_b'] 

	else: 
		#print x 
		snpa=x[2]; snpb=x[5] 
		if snpa in map_data: x+=[map_data[snpa]]
		else: x+=['na']
		if snpb in map_data: x+=[map_data[snpb]]
		else: x+=['na']

	print "\t".join(x) 	
		
