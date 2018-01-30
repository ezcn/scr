d_snps={}
for line in open ("/home/enza/oogaprotocol/PEPERONE/9_domesticationpepper/bedintersections/all.intersection_domreg_vcfkim.counts", "r"): 
	x=line.rstrip().split()
	xkey=x[0].split(".")[0]
	d_snps[xkey]=x[1]

d_coordkim={}
for line in open ("/home/enza/oogaprotocol/PEPERONE/9_domesticationpepper/dataset17_qin2014.blast.kim2014.besthit.bed", "r" ) : 
	y=line.split() 
	ykey=y[3].split(".")[0]
	d_coordkim[ykey]=[y[0], y[1], y[2]]



#print d_snps
#print d_coordkim

#"""
count =0
for line in open ("/home/enza/oogaprotocol/PEPERONE/9_domesticationpepper/dataset17_qin2014.txt", "r") : 
	z=line.rstrip().split("\t") 
	if count ==0 : 
		res=["chrkim", "startkim", "endkim", "snpGBS"]+z 
		print "\t".join(res )
		count +=1 
	else: 
		if z[0] in d_snps: 
			res=d_coordkim[z[0]]+[d_snps[z[0]]] +z 		
			print "\t".join(  map(str, res)  ) 
		else: 
			res=["na", "na", "na", "na"] +z
			print "\t".join(  map(str, res)  )
#"""
