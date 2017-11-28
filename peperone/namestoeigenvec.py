dname={}

for line in open("/home/enza/oogaprotocol/PEPERONE/pop2/list.complete", "r") : 
	x=line.split()
	dname[x[0]]=x

for line in open("/home/enza/oogaprotocol/PEPERONE/1_pca/gw.eigenvec", "r") : 
	y=line.split() 
	print "\t".join( dname[y[0]]+y) 

