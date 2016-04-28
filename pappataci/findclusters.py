import sys
import re 

if len(sys.argv) < 2:
    print "\nusage", sys.argv[0], "input_file", "min_hit", "cl_len", "step"
    print "~~~~~~~~~~" 
    print "#  find cluster of regex in  sequences exploring sequence in sliding windows\n" 
    print "input_file: see test.in" 
    print "minhit: minimum number of regex in cluster" 
    print "cl_len: length of the cluster in bp "
    print "step: step of the sliding window\n"

    exit(1)


infile=sys.argv[1] # input file 
minhits=int(sys.argv[2]) #numero minimo di regex in un cluster 
clusterlen=int(sys.argv[3]) #lunghezza del cluster 
step=int(sys.argv[4]) # passo sliding windows 


#~~~~~~~~~~~~~~~

seq_dz={}

count=0
for line in open(infile, "r") :
	if  re.search("#--" , line) : count +=1   
	x=line.rstrip().split()
	if re.search("Sequence:" , line): 
		seqname=x[2]; sfrom=int(x[4]);  sto=int(x[6])
		seq_dz[seqname]=[[sfrom ,sto] ]
	if re.search("regex", line): 
		rfrom=int(x[0]); rto=int(x[1]); rseq=x[4]
		seq_dz[seqname].append([rfrom, rto, rseq] ) 

"""
#print seq_dz

title=["seqname", "seqstart", "seqend","hitcount", "regstart", "regend", "regseq"]
print "\t".join(title) 

for ss in seq_dz:
	for rr in seq_dz[ss][1:]: 
		res=[ss] + seq_dz[ss][0] +[len(seq_dz[ss]) -1 ] + rr 
		print "\t".join(map(str, res) )
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cluster_dz={}
for ss in seq_dz:
	cluster_dz[ss]=[]
	seqlen = seq_dz[ss][0][1] - seq_dz[ss][0][0]	
	if seqlen < clusterlen: pass # seq piu piccola del cluster 
	elif  len(seq_dz[ss]) -1  < minhits : pass # hits minori del minimo richiesto 
	else: 
		for y in range(0, seqlen, step):
			start=y ; end=y+clusterlen; clulist=[]
			for reg in seq_dz[ss][1:]:
				if reg[0]>=start and reg[1]<=end: 
					clulist.append(seq_dz[ss].index(reg)) 
			if len(clulist) >= minhits: 
				cluster_dz[ss].append(clulist) 
#print cluster_dz

for seq in cluster_dz: 
	if len(cluster_dz[seq]) >1 : # se almeno un cluster
		# lista unica 
		templist=[]
		for i in cluster_dz[seq]: 
			if not i in templist: templist.append(i) 

		# recupera dati 
		countcl=1
		for cluster in templist:
			res=["%s_cl%s" %(seq , countcl), len(cluster)]
			for elementindex in cluster:  
				res.append( seq_dz[seq][elementindex]) 
			countcl+=1 
			print res   


		
