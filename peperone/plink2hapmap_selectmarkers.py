import csv
import sys 

mysiteslist=sys.argv[1] 
mymap=sys.argv[2]
mytransposedped=sys.argv[3]

#~~~~~~  read marker list 
sitestoinclude=[]
for iline in open(mysiteslist) : 
	y=iline.split() 
	sitename=y[1]
	sitestoinclude.append(sitename ) 		

#~~~~~   read map file 
linestokeep=[]		
d_sites={}
countline=1
for mline in open (mymap) : 
	m=mline.split() 
	tempsitename=m[1]
	if tempsitename in sitestoinclude: 
		linestokeep.append(countline)
		d_sites[countline]=m 
	countline+=1 

#~~~~~~  Read transposed ped file and write a hapmap output 
counttrline=0 
for tline in  open(mytransposedped):
	#print tline 
	t=tline.rstrip().split()
	if counttrline==0 : 
		title=['rs', 'alleles', 'chrom', 'pos', 'strand', 'assembly', 'center', 'protLSID', 'assayLSID', 'panelLSID','QCcode']
		title+=t 
		print " ".join(map(str, title ))  

	else:  
		if counttrline in linestokeep: 
			templine=[]
			[templine.append("na") for f in range(11)] 
			templine[0]=d_sites[counttrline][1]
 			templine[2]=d_sites[counttrline][0]
			templine[3]=d_sites[counttrline][3]
			templine+=t 
			print " ".join(map(str, templine ))  
	counttrline+=1 



