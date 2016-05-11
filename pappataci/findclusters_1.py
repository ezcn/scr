import sys
import re 

#TEST:  python findclusters_1.py  test.dreg 3 3000 50

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

#print seq_dz

"""
title=["seqname", "seqstart", "seqend","hitcount", "regstart", "regend", "regseq"]
print "\t".join(title) 

for ss in seq_dz:
	for rr in seq_dz[ss][1:]: 
		res=[ss] + seq_dz[ss][0] +[len(seq_dz[ss]) -1 ] + rr 
		print "\t".join(map(str, res) )
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find clusters according to input parameters and store info in cluster_dz 
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# print header 
head=["seq", "cluster", "nb_elements" , "cluster_lenght", "cluster_start", "cluster_end", "elements"]
print "\t".join(head) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#retrieve sequence info for every transcript (seq)
for seq in cluster_dz: 
	if len(cluster_dz[seq]) >1 : # transcripts with at least one cluster
		# make a list of unique clusters 
		uniquelist=[]
		for i in cluster_dz[seq]: 
			if not i in uniquelist: uniquelist.append(i) 
		#print "*****", uniquelist 

		# get data for clusters and store in clusterdata 
		countcl=1; clusterdata={}
		for cluster in uniquelist:
			templist=[]
			for elementindex in cluster:
				templist.append( seq_dz[seq][elementindex])
			cluster_start=templist[0][0]
			cluster_end=templist[-1][1]
			cluster_length=cluster_end-cluster_start
			tempkey=("%s_cl_%s" %(seq , countcl), len(cluster) , cluster_length, cluster_start, cluster_end )

			clusterdata[tempkey]=templist
			countcl+=1
		#for i in clusterdata: print i, clusterdata[i]
		#print "************************************" 

		# find nested clusters 	
		numberofclusters=len(clusterdata.keys()) 
		lunghi=[]
		while  numberofclusters >0 : 

			# identifica cluster piu lungo (to fix: se ci sono due cluster di uguale lunghezza non discrimina) 
                	lengths=[i[2] for i in clusterdata.keys()]
                	longestcluster=max(lengths)

			#first iteration to find longest transcript 
			for cl in clusterdata.keys(): 
				#print "''''''''''''''''", cl   
				if cl[2]==longestcluster: 
					if not (cl, clusterdata[cl]) in lunghi: 
						lunghi.append((cl, clusterdata[cl]));
						#print "longesttoberemoved" # , "lllllll" 
						longeststart=clusterdata[cl][0][0]
						longestend=clusterdata[cl][-1][1]	
						del clusterdata[cl]
						numberofclusters-=1 	

			# second iteration to find nested clusters 
			for  clcheck in clusterdata.keys(): 
				#print "++++++" , clcheck 	
				if clcheck[2]<longestcluster: 
					checkstart=clusterdata[clcheck][0][0]
					checkend=clusterdata[clcheck][-1][1]
					if checkstart >= longeststart and  checkend <= longestend: 
		   				del clusterdata[clcheck]							
						numberofclusters-=1
	
		#format output 
		for r in lunghi:
			resline=[]; clusterstring=""
			resline.append(seq)  
			for ii in r[0]: resline.append(ii)	
			for cc in map(str, r[1]): 
				for ee in cc: 
					clusterstring+=ee  #lunghi  
			resline.append( clusterstring)
			print "\t".join(map(str, resline))

				
