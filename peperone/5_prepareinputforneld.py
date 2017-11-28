import re 

d_pos={}
for c in range (1 , 13) : 
	infile="/home/enza/oogaprotocol/PEPERONE/geneticmap/%s.phys.gen.map" % (c) 
	#print infile 
	with open (infile, "r") as f: 
		next(f) 
		for line in f : 
			if not re.search('NaN', line):  
				x=line.split() 
				snp=x[2]; phy=float(x[3]); gen=float(x[5].rstrip(")") )   #gen
				d_pos[snp]=[phy, gen ] 

title=['chra', 'bpa', 'snpa', 'chrb', 'bpb', 'snpb', 'r2','type', 'bpa-bpb', 'phy_ab' , 'gen_ab']
print "\t".join(title) 

with open ("/home/enza/oogaprotocol/PEPERONE/3_ld/all.ldout") as g: 
	next (g) 
	for gline in g: 
		y=gline.split() 
		snp1=y[2]; snp2=y[5]; pos1=int(y[1]); pos2=int(y[4]) 
	
		if snp1 in d_pos and snp2 in d_pos: 	
			res=y + [abs(pos1-pos2)] + [abs(x - y)  for x, y in zip(d_pos[snp1], d_pos[snp2])]
			print "\t".join(map(str, res ) )  
