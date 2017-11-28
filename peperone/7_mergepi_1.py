import sys, re 


print "\t".join(['chr', 'pos' , 'pi', 'type' ]) 
types=['wild', 'domesticated']
for  c in range (1, 13): 
	for t in types: 
		pathtofile="/home/enza/oogaprotocol/PEPERONE/7_pi/%s.pi.%s.sites.pi" % (c,  t)  
		with open (pathtofile, "r" ) as f: 
			next(f) 
			for line in f:
				x=line.split()
				res=x+[t]
				print "\t".join( map(str, res) )  
		 

						
