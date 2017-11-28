import os 
import sys 


cla={}
with open("/home/enza/oogaprotocol/PEPERONE/pop2/list.clades") as c:
	clades = c.read().splitlines()	 
	for cc in clades: cla[cc]=[]

spe={}
with open("/home/enza/oogaprotocol/PEPERONE/pop2/list.species") as s:
	species = s.read().splitlines()
	for ss in species: spe[ss]=[]

typ={}
with open("/home/enza/oogaprotocol/PEPERONE/pop2/list.types") as t:
	types = t.read().splitlines()
	for tt in types: typ[tt]=[]

with open("/home/enza/oogaprotocol/PEPERONE/pop2/list.complete") as l: 
	for line in l: 
		x=line.split(); sp=x[1]; cl=x[4]; ty=x[5]; missingness=x[6]	
		if missingness=="ok": 
			cla[cl].append(x)
			spe[sp].append(x) 
			typ[ty].append(x) 

"""
print cla.keys() 
print typ.keys() 
print spe.keys()
"""
for cname in cla: 
	out=open("/home/enza/oogaprotocol/PEPERONE/pop2/pop.cla.%s" %(cname), "w") 
	sys.stdout=out 
	for item in cla[cname]: 
		print "\t".join(item ) 


for sname in spe:
        out=open("/home/enza/oogaprotocol/PEPERONE/pop2/pop.spe.%s" %(sname), "w")
        sys.stdout=out
        for item in spe[sname]:
                print "\t".join(item )

for tname in typ:
        out=open("/home/enza/oogaprotocol/PEPERONE/pop2/pop.typ.%s" %(tname), "w")
        sys.stdout=out
        for item in typ[tname]:
                print "\t".join(item )
