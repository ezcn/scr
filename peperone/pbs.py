import math, sys 

fstfileSN=sys.argv[1]
fstfileSO=sys.argv[2]
fstfileNO=sys.argv[3] 



def PBS( fstSN, fstSO, fstNO):
	"population branching statistics "
	tSN=-math.log(1.00000001-fstSN, 10)
 	tSO=-math.log(1.00000001-fstSO, 10)
	tNO=-math.log(1.00000001-fstNO, 10) 
	PBS=(tSN + tSO -tNO)/2 
	return PBS 

def range01 (myvalue, mylist) :  
	lowest=min(mylist) ; highest=max(mylist) 
	myrangevalue= (myvalue-lowest)/(highest-lowest)	
	return myrangevalue 


d_fst={}

filelist=[fstfileSN, fstfileSO, fstfileNO]
 
for ffile in filelist: 
	with open(ffile) as f:
		next(f)
		for line in f:
			x=line.split()
			mychr=x[0];  mypos=x[1]; myfst=x[2]; fst=0 
			if myfst=="-nan": continue 
			elif float(myfst)<0 : continue  #fst=0 
			else : fst=float(myfst) 
			if not (mychr, mypos) in d_fst: d_fst[(mychr, mypos)]=[]
			d_fst[(mychr, mypos)].append(fst )


#print d_fst 
title =['chr', 'pos', 'fstSN', 'fstSO', 'fstNO', 'PBS']	
print "\t".join(title )		

#d_pbs={}
#listpbs=[]

for i in d_fst: 
	if len(d_fst[i])==3: 
		#if not (0 in d_fst[i]  or 1 in d_fst[i]) : 
		res=[i[0], i[1]]
		mypbs= PBS(d_fst[i][0], d_fst[i][1], d_fst[i][2])
		if mypbs > 0: 
			res+=d_fst[i]
			res.append(mypbs) 
			print "\t".join(map(str, res ) ) 
			#listpbs.append(mypbs) 
			#d_pbs[i]=mypbs 
			#print i[0], i[1], mypbs  
		#else: print d_fst[i]


#title =['chr', 'pos', 'fstSN', 'fstSO', 'fstNO', 'PBS', 'PBS01']	
#print "\t".join(title )		
#for p in d_pbs: 
#	norpbs=range01(d_pbs[p], listpbs  )	
#	res=[p[0], p[1]]+d_fst[p]
#	res.append(d_pbs[p])	
#	res.append(norpbs) 
#	print "\t".join(map(str, res ) )



