import sys, re 


print "\t".join(['chr', 'pos' , 'pi_wild', 'pi_domesticated' ]) 
types=['wild', 'domesticated']
for  c in range (1, 13): 
	data_pi={} ; pos_ord=[] 
	for t in types: 
		pathtofile="/home/enza/oogaprotocol/PEPERONE/7_pi/%s.pi.%s.sites.pi" % (c,  t)  
		with open (pathtofile, "r" ) as f: 
			next(f) 
			for line in f:
				#print line 
				x=line.split() 
				chrom=x[0]; pos=x[1]; pi=x[2] 
				pos_ord.append(pos) 
				if not pos in data_pi: data_pi[pos]={}
				data_pi[pos][t]= pi
	for pp in pos_ord: 
		res=[c, pp]
		for tt in types: res.append(data_pi[pp][tt]) 
		print "\t".join( map(str, res) )  
		 

						
