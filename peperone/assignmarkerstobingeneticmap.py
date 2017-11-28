import pandas as pd 
import re 
import sys 

geneticmapfile=sys.argv[1]
physicalmapfile=sys.argv[2]
###  works one chr at time

##### define bin and bin names
bins=[0]
binnames=[]
with open(geneticmapfile, "r" ) as b:  
	next(b) 
	for bline in b: 
		y=bline.split() 
		bins.append(float(y[2])); binnames.append((y[1], y[3]) )  

#del binnames[-1] 
#print binnames 


######  input file  with position to assing 
d_index={0:"chr", 1: "snp", 3:"phys"}
d_map={"chr": [], "snp":[], "phys":[] }
colname=["chr", "snp", "phys"]

with open(physicalmapfile, "r" ) as f:  
	for line in f: 
		x=line.split() 
		d_map["chr"].append(x[0]) 
		d_map["snp"].append(x[1]) 
		d_map["phys"].append(float(x[3]))  
#		for f in x :i
#			if x.index(f) != 2:  
#				d_map[ d_index[x.index(f)] ].append(f)  	

#print d_map 
#print colname 

df = pd.DataFrame(d_map, columns =colname) 
df['categories'] = pd.cut(df['phys'], bins, labels=binnames)

#print df 
with pd.option_context('display.max_rows', None, 'display.max_columns', 5):
    print(df)
