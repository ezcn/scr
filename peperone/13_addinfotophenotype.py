
d_infocomplete={}

for line in open ("/home/enza/oogaprotocol/PEPERONE/pop4/list.complete.pop4"): 
	x=line.rstrip().split('\t') 
	d_infocomplete[x[0]]=x


d_q={}
for yline in open ("/home/enza/oogaprotocol/PEPERONE/pop4/list.species.annuum.q1.q7"): 
	y=yline.rstrip().split('\t')
	#print y 
	d_q[y[0]]=y[1]

#print d_q 

countme=0
for zline in open ("/home/enza/oogaprotocol/PEPERONE/13_pheno/pheno.data"): 
	z=zline.rstrip().split('\t') 
	if countme==0: 
		countme+=1 
		title=z+['ID','species', 'country', 'continent', 'clade', 'dom_wild', 'ok_missing', 'q' ]
		print "\t".join(title) 
	else: 
		myid=z[1]; res=z
		if myid in d_infocomplete: res+=d_infocomplete[myid]
		else: res+=["na" for i in range(7)]
		if myid in d_q: res.append(d_q[myid]) ; #print res 
		else:  res+=['q0']

		print "\t".join(map(str, res ))



	
