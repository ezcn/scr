
d_chr={}
for line in open('/home/enza/oogaprotocol/PEPERONE/0_GBS/chrlength.tsv'): 
	x=line.split() 
	mychr=x[1]; chrlen=x[3]
	d_chr[mychr]=chrlen 

#print d_chr 
count=0 
for fline in open('/home/enza/oogaprotocol/PEPERONE/0_GBS/pepper.master_tags.sorted.gbs_captures.cov.snps.dist.tsv'): 
		
	y=fline.rstrip().split('\t') 

	if count==0: 
		y.append('chr_size')
		count+=1 
	
	else: y.append(d_chr[y[0]])

	print  "\t".join(map(str, y  ))
	 
