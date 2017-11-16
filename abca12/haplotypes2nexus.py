import random, re, sys 

haplotypefile=sys.argv[1]

d_pop={}
#for line in open ("datafileshaplotypes/all.complete.samples", "r" ) : 
for line in open ("datafileshaplotypes_enza/samples.info.all.list", "r" ) : 
	x=line.rstrip().split(); idsample = x[0] 
	d_pop[idsample]=x[2]


#print d_names


d_da={0:{0:"D", 1:"A"}, 1:{0:"A", 1:"D"}, 2:{0:"A", 1:"D"}, 3:{0:"D", 1:"A"}, 4:{0:"A", 1:"D"}, 5:{0:"A", 1:"D"}, 6:{0:"A", 1:"D"}, 7:{0:"D", 1:"A"}, 8:{0:"D", 1:"A"}, 9:{0:"D", 1:"A"}} 


l_pop=[]; d_fasta={}; d_seqcounter={}
seqcounter=1
for line in open (haplotypefile, "r"): 
	y=line.rstrip().split("\t") 
	#myid=y[0]; mypop=d_pop[myid]; myseq=y[1].replace(" ", "").replace(".", "?") 
	myid=y[0][:-2]; mypop=d_pop[myid]; myseq=y[1]
	mydaseq=""; myseqpos=0
	for i in myseq: 
		mydaseq+=d_da[myseqpos][int(i)]
		myseqpos+=1
	l_pop.append(mypop) 
	if not myseq in d_fasta: 
	#if not mydaseq in d_fasta: 
		d_fasta[myseq]=(seqcounter, [mypop])
		#d_fasta[mydaseq]=(seqcounter, [mypop])			
		d_seqcounter[seqcounter]=(myseq, mydaseq) 
		#d_seqcounter[seqcounter]=mydaseq
		seqcounter+=1 
	else: 
		d_fasta[myseq][1].append(mypop) 
		#d_fasta[mydaseq][1].append(mypop)

#print d_fasta# ["T A T G T T T A A G A"]	



###################################
minNbSeq=10 # minimun frequency of sequence 
ntaxa=0
highdsitepos=9
for i in d_fasta: 
	 if len(d_fasta[i][1])>= minNbSeq: ntaxa+=1

#print "#NEXUS\nBEGIN TAXA;\nDIMENSIONS NTAX=%s;\nTAXLABELS" %( ntaxa)
#for i in d_fasta.keys(): 
#for i in range(1, seqcounter):
#	if len(d_fasta[d_seqcounter[i]][1])>minNbSeq:  #check if frequency of seq is greater that minNbSeq 
#		print "seq_%s " % (i )#(d_fasta[i][0]) 
#print ";\nEND;\n\n" 

seqlen=len(random.choice(list(d_fasta))) 
#print "BEGIN CHARACTERS;\nDIMENSIONS NCHAR=%s;\nFORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. ;\nMATRIX" % ( seqlen) 
#for j in d_fasta.keys(): 
#for j in range(1, seqcounter): 
#	if len(d_fasta[d_seqcounter[j]][1])>minNbSeq: #check if frequency of seq is greater that minNbSeq 
#		print "seq_%s %s" % (j , d_seqcounter[j]) #(d_fasta[j][0] , j) 
#print ";END;\n\n"

print "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=%s;" %( ntaxa)    #\nTAXLABELS" %( ntaxa) 
print "DIMENSIONS NCHAR=%s;\nFORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. ;\nMATRIX" % ( seqlen)
for j in range(1, seqcounter): 
	if len(d_fasta[d_seqcounter[j][0]][1])>minNbSeq: #check if frequency of seq is greater that minNbSeqi
		highdallele=d_seqcounter[j][0][highdsitepos]
		myletter="A" if highdallele=="1"  else "G" 
		print "seq_%s_%s %s" % (d_seqcounter[j][1], j , d_seqcounter[j][0]) #(d_fasta[j][0] , j) 
print ";\nEND;\n\n"


print "BEGIN TRAITS;\n Dimensions NTRAITS=%s;\n Format labels=yes missing=? separator=Comma;\nTraitLabels" % (len(set(l_pop)) ), 
for pop in set(l_pop) : print pop,  

print ";\nMatrix" 
#for seq in d_fasta: 
for s in d_seqcounter: 
	if len(d_fasta[d_seqcounter[s][0]][1])>minNbSeq: #check if frequency of seq is greater that minNbSeq 
		highdallele=d_seqcounter[s][0][highdsitepos]
		myletter="A" if highdallele=="1"  else "G" 
		print "seq_%s_%s" % (d_seqcounter[s][1] , s), #(d_fasta[seq][0]),  
		tempocc=""
		for p in set(l_pop):
			myoccurrence= d_fasta[d_seqcounter[s][0]][1].count(p)
			tempocc+=str(myoccurrence)
			tempocc+=','  
		striptemocc=tempocc.rstrip(',') 
		print striptemocc  
	#print "\r"	

print ";\nEND;"



#matching characters same as missing """" 
