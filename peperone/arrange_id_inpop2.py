
with open('/home/enza/oogaprotocol/PEPERONE/pop2/list.missing.ids') as m:
    missing = m.read().splitlines()
#print missing 

d_clades={}
with open("/home/enza/oogaprotocol/PEPERONE/pop2/list.clade") as f:
    next(f)
    for line in f:
	x=line.split()	
	clade =x[0] ; species=x[1] ; tipo=x[3]
	d_clades[species]=[clade, tipo]
	d_clades["failed"]=["failed", "failed"]

with open ("/home/enza/oogaprotocol/PEPERONE/pop2/list.species.geographic.origin") as g: 
	next(g) 
	for line in g: 	
		y=line.split()
		ind=y[3]; spec=y[2]; 
		datamiss="missing" if ind in missing else "ok"
		res=[ind, spec] +y[0:2]+ d_clades[spec] 
		res.append(datamiss) 
		print "\t".join(res)  

			
 
