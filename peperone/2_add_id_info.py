import sys

k=int(sys.argv[1]) 
infile=sys.argv[2] 

id_data={} 
for line in open ("/home/enza/oogaprotocol/PEPERONE/pop2/list.complete"):  
	x=line.rstrip().split() 
	id_data[x[0]]=x

title=['id', 'clade', 'nat', 'con', 'species', 'type', 'miss', 'id2']
for i in range (1, k+1) : 
	title.append("c%s"  % (i)) 

print  "\t".join(title )


for ine in open (infile, "r") : 
	y=ine.rstrip().split()
	res=id_data[y[0] ] +y 
	print  "\t".join(res) 
	   
