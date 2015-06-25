#!/bin/python
import numpy as np
import re 
import sys 
   

pop=sys.argv[1]
file1=sys.argv[2] # relat.v31 lista delle coppie 
file2=sys.argv[3] # YOB
file3=sys.argv[4] # kinship 
file4=sys.argv[5]  # q statistis (pedigree informativeness 
file5=sys.argv[6] # inbreeding 

f1=open(file1, "r")
f2=open(file2, "r")
f3=open(file3, "r")
f4=open(file4, "r") 
f5=open(file5, "r") 

yearbins=range(1525, 1965, 5)  
kinbin=[0.031250000, 0.007812500, 0.001953125, 0.000488281, 0.000122070, 0.000030518, 0.000]
############################################################
def assignbin(mylist, myvalue): 
	mybin=len(mylist)+1 ; mylistv=0
	for v in mylist:  
		if myvalue>=v : mybin-=1; mylistv=v
	return (mybin, mylistv) 

###########################################################
mycouple_list=[] ; mydic_nfigli={}; mydic_year={};  mydic_kin={}; mydic_idchi={}; mydic_parents={};  mydic_chir={} ; mydic_grchi={}; mydic_q4={} ; mydic_q8={} 
av_bin_chi={};  av_bin_chir={}; av_bin_grchi={}; inbreeding_dic={}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


count_yobisna=0
for line2 in f2 : #anno dell'individuo 
	if not re.search("PID", line2) : 
       		s=line2.split()
		id=s[0]; ystart=int(s[2]) ; yend=int(s[3]) 
		if ystart==yend: 
			yob=ystart 
		elif abs(ystart-yend)<=10 : 
			yob=(ystart+yend)/2 
		else: yob="na" ; count_yobisna+=1 	
		mydic_year[id]=yob
#print mydic_year 
for ll in f4: 
	q=ll.rstrip().split()
	mydic_q4[q[0]]=float(q[1]) 
	mydic_q8[q[0]]=float(q[3]) 
	
for line in f1 :   #coppie e figli  
	t=line.split(",") 
	couple=(t[0], t[1]) 
	mycouple_list.append(couple)
	nfigli=t[3].count('<><>')+1
	mydic_nfigli[couple]=nfigli 

	if not t[0] in mydic_nfigli: 
		mydic_nfigli[t[0]]=nfigli 
	else: mydic_nfigli[t[0]]+=nfigli

	if not t[1] in mydic_nfigli: 
		mydic_nfigli[t[1]]=nfigli
	else: mydic_nfigli[t[1]]+=nfigli

	u=line.rstrip().split()[1:][0::2] 
	mydic_idchi[couple]=[xu.rstrip().lstrip() for xu in u]

for ccp in mycouple_list: 
	fa=ccp[0]; mo=ccp[1]
	yofa=mydic_year[fa]; yomo=mydic_year[mo]
	if not (yofa=='na' or yomo =="na") :
		yoc=(yofa+yomo)/2	
	else: yoc="na"
	mydic_year[ccp]=yoc

for cp in mycouple_list: ### nb childern that reproduce  and grandchilden
	count_chir=0; count_grchi=0
	for ch in mydic_idchi[cp]:
		if ch in mydic_nfigli:
			if int(mydic_nfigli[ch])>0 : count_chir+=1; count_grchi+=int(mydic_nfigli[ch])
	mydic_chir[cp]=count_chir ; mydic_grchi[cp]=count_grchi

list_chi_parentskin_is_na=[]
for line3 in f3 : #kinship e kinship bins 
        z=line3.split()
      	zz=z[0].split("_") 
	couple3=(zz[0], zz[1])
	kin=z[1]
	if not kin=='na': 
		bin=assignbin(kinbin, float(kin)) 
        	mydic_kin[couple3]=( kin, bin[0]) 
	else: 
		mydic_kin[couple3]=("na", "na") 
		for i in mydic_idchi[couple3]: 
			list_chi_parentskin_is_na.append(i) 

for line5 in f5: #inbreeding 
	n=line5.rstrip().split()
	nind=n[0]; ninbrd=n[2]
	
	if not nind in list_chi_parentskin_is_na: 
		inbreeding_dic[nind]=ninbrd
	else: inbreeding_dic[nind]="na"

#~~~~~~~~~ STANDARDIZATION ~~~~~~~~~~~~~~~~~
for b in yearbins: 
	av_bin_chi[b]=[];  av_bin_chir[b]=[]; av_bin_grchi[b]=[]

#print av_bin_chi

c_kin=0.0000001; c_chi=1  
for cc in mycouple_list:
	#print "#############", cc
 	yoc=mydic_year[cc]
	if not  yoc=="na" : 
		bb=assignbin(yearbins, yoc)
		#print yoc, bb[1]  	
		nchi=mydic_nfigli[cc];   chir=mydic_chir[cc];   grchi=mydic_grchi[cc]
		#print nchi,  chir, grchi 

		p_nchi=np.log(nchi+c_chi); p_chir=np.log(chir+c_chi);   p_grchi=np.log(grchi+c_chi)
			
		av_bin_chi[bb[1]].append(p_nchi) 
		av_bin_chir[bb[1]].append(p_chir)
		av_bin_grchi[bb[1]].append(p_grchi)

#print "#individuals with no yob:" , count_yobisna

mytitle=["pop", "couple", "kin", "bin" , "yfa", "ymo", "yoc", "nchi","chir", "grchi", "snchi", "schir", "sgrchi", "iil_av" , "iil_sd", "moatf", "faatf" , "q4", "q8", "inb_fa", "inb_mo"]
print "\t".join(mytitle)

for c in mycouple_list: 
	fa=c[0]; mo=c[1]
	yofa=mydic_year[fa]; yomo=mydic_year[mo]

	### nchi , chir, grchi, yoc , kin 
	nchi= int(mydic_nfigli[c]) ; chir=mydic_chir[c];   grchi=mydic_grchi[c] 
	yoc=mydic_year[c] 

	if c in mydic_kin:  
		if not mydic_kin[c][0]=="na": 
			ckin=float(mydic_kin[c][0])+c_kin ;  cbin=mydic_kin[c][1]
		else: ckin="na" ; cbin="na"
	else: ckin="na"	; cbin="na"		
	
	### standardized variables 
	if not  yoc=="na" :
		bb=assignbin(yearbins, yoc)
		p_nchi=np.log(nchi+c_chi); p_chir=np.log(chir+c_chi);   p_grchi=np.log(grchi+c_chi)

		snchi=( p_nchi - np.mean(av_bin_chi[bb[1]]) ) /float( np.std(av_bin_chi[bb[1]]) )  
		snchir=( p_chir - np.mean(av_bin_chir[bb[1]]) ) /float( np.std(av_bin_chir[bb[1]]) ) 
		sngrchi=( p_grchi - np.mean(av_bin_grchi[bb[1]]) ) /float( np.std(av_bin_grchi[bb[1]]) ) 

	else: snchi=snchir=sngrchi="na"	

	###  Interbirth Intervals Length, age at first child 
	yearchi=[]
	#print c, mydic_idchi[c]
	for chichi in mydic_idchi[c]: 
		if not mydic_year[chichi]=="na": 
			yearchi.append(int(mydic_year[chichi])) 
	if len(yearchi)>1: 
		yearchi.sort()
		#print yearchi, yomo, yofa
		myinterlist=[yearchi[x+1]-yearchi[x] for  x in range(len(yearchi)-1) ]
		iil_av=np.mean(myinterlist)
		#iil_av=sum( myinterlist)/float(len(myinterlist)) 
		iil_sd=np.std(myinterlist)

	else : iil_av=iil_sd="na"

	if len(yearchi)>0:
		if yomo !="na": moatf=yearchi[0]-int(yomo)
		else: moatf="na" 	
		if yofa !="na": faatf=yearchi[0]-int(yofa) 
		else:  faatf="na" 	
	else : moatf=faatf="na" 

	q4=q8="na"
	if  (fa in mydic_q4 and  mo in mydic_q4 ) : q4=mydic_q4[fa]*mydic_q4[mo]
	if  (fa in mydic_q8 and  mo in mydic_q8 ) : q8=mydic_q8[fa]*mydic_q8[mo]

	myoverlap=[ pop, "%s_%s" %(fa, mo), ckin, cbin, yofa, yomo, yoc , nchi, chir, grchi, snchi, snchir, sngrchi, iil_av , iil_sd, moatf, faatf , q4, q8, inbreeding_dic[fa], inbreeding_dic[mo] ]
	#print myoverlap 
	print "\t".join(str(x) for x in myoverlap)


