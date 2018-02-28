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

title=["CHR_A", "BP_A", "SNP_A", "CHR_B","BP_B", "SNP_B", "R2", "species", "clade", "type", "BP_A-BP_B", "phy_ab", "gen_ab"]
print "\t".join(title) 

with open ("/home/enza/oogaprotocol/PEPERONE/3_ld_species_q/all.species.ld") as g: 
	next (g) 
	for gline in g: 
		y=gline.split() 
		snp1=y[2]; snp2=y[5]; pos1=int(y[1]); pos2=int(y[4]) 
	
		if snp1 in d_pos and snp2 in d_pos: 	
			res=y + [abs(pos1-pos2)] + [abs(x - y)  for x, y in zip(d_pos[snp1], d_pos[snp2])]
			print "\t".join(map(str, res ) )  
