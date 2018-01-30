import random, re, sys 

haplotypefile=sys.argv[1]
"""
A418:C5W6FACXX:2:250422395-1    000000
A98:C5W6FACXX:2:250422392-0     000000
A98:C5W6FACXX:2:250422392-1     000000
AjiCristal:C5W6FACXX:2:250422367-0      000000
AjiCristal:C5W6FACXX:2:250422367-1      000000
AjiGinger:C5W6FACXX:3:250422590-0       000000
AjiGinger:C5W6FACXX:3:250422590-1       000000

"""

d_pop={}; d_wild_dom={}
l_latitude=[]
for line in open ("/home/enza/oogaprotocol/PEPERONE/pop4/list.complete.pop4", "r" ) : 
	x=line.rstrip().split(); idsample = x[0] ; #latitudeinterval=x[5]
	mypop=x[1]; my_wild_dom=x[5]
	d_pop[idsample]=x[1]
		
	if not mypop in d_wild_dom: d_wild_dom[mypop]=my_wild_dom
#################################


l_pop=[]; d_fasta={}; d_seqcounter={}
seqcounter=1
for line in open (haplotypefile, "r"): 
	y=line.rstrip().split("\t") 
	myid=y[0].split("-")[0]
	if myid in d_pop: 
		mypop=d_pop[myid]; myseq=y[1]
		#print myid, mypop, myseq 
		l_pop.append(mypop) 
		if not myseq in d_fasta: 
			d_fasta[myseq]=(seqcounter, [mypop])
			d_seqcounter[myseq]=seqcounter		
			seqcounter+=1 
		else: 
			d_fasta[myseq][1].append(mypop) 

#print d_fasta # ["T A T G T T T A A G A"]	
#print d_wild_dom

#"""
###################################
minNbSeq=2 # minimun frequency of sequence 

for sequence in d_fasta: 
	#print '################'
	if len( d_fasta[sequence][1] )>minNbSeq: #check if frequency of seq is greater that minNbSeqi
		temppop=set(d_fasta[sequence][1])
		for p in temppop:
			#print p, d_wild_dom[p] 
			popoccurrence= d_fasta[sequence][1].count(p)	
			#print popoccurrence
			print ">s%s_%s_%s_%s" % (d_seqcounter[sequence], p, d_wild_dom[p], popoccurrence ) #(d_fasta[j][0] , j) 
			print sequence
"""
print ";\nEND;\n\n"


print "BEGIN TRAITS;\n Dimensions NTRAITS=%s;\n Format labels=yes missing=? separator=Comma;\nTraitLabels" %  (len(set(l_pop)) ), #len(set(l_latitude))#(len(set(l_pop)) ), 


for pop in set(l_pop) : print pop,  
#for ll in set(l_latitude): print ll 

print ";\nMatrix" 
#for seq in d_fasta: 
for s in d_seqcounter: 
	if len(d_fasta[d_seqcounter[s]][1])>minNbSeq: #check if frequency of seq is greater that minNbSeq 
		print "seq_%s" % ( s), #(d_fasta[seq][0]),  
		tempocc=""
		for p in set(l_pop):
		#for p in set(l_latitude): 
			myoccurrence= d_fasta[d_seqcounter[s]][1].count(p)
			tempocc+=str(myoccurrence)
			tempocc+=','  
		striptemocc=tempocc.rstrip(',') 
		print striptemocc  
	#print "\r"	

print ";\nEND;"



""" 
