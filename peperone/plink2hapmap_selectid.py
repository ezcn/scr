import csv
import sys 

myped=sys.argv[1]
myidlist=sys.argv[2] 
myout=sys.argv[3]


#~~~~~~  read id list 
indvtoinclude=[]
for iline in open(myidlist) : 
	y=iline.split() 
	indvtoinclude.append(y[0]) 		



#~~~~~~  Read ped file and write a transposed 
mylistofidgenotypes=[] 
for line in  open(myped):
 	x=line.rstrip().split()
	templist=[]; indvID=x[0]
	if indvID in indvtoinclude: 
		templist.append(indvID) 
		for i in range(6, len(x), 2) : 
			genotype=''.join(map(str, [x[i], x[i+1]] ) ) 
			templist.append(genotype )
		mylistofidgenotypes.append(templist)

#print mylistofidgenotypes

transposed = zip(*mylistofidgenotypes)
with open(myout , 'w') as fout:
	w = csv.writer(fout, delimiter='\t')
	w.writerows(transposed)

#

