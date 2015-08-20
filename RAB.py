import sys,  random , numpy , scipy.stats

myinputderivedallelecounts=sys.argv[1]
mypop1haploidsize=int(sys.argv[2]) 
mypop2haploidsize=int(sys.argv[3])
mylistofsites=sys.argv[4]
mygroups=int(sys.argv[5]) 


"""
ref  Do R, Balick D, Li H, Adzhubei I, Sunyaev S, Reich D. No evidence that
selection has been less effective at removing deleterious mutations in Europeans 
than in Africans. Nat Genet. 2015 Feb;47(2):126-31. doi: 10.1038/ng.3186. Epub
2015 Jan 12. PubMed PMID: 25581429; PubMed Central PMCID: PMC4310772.


1) myinputderivedallelecounts -> chromosome position counts_in_pop_1 counts_in_pop_2 
four  columns  no header space or tab delimited

2) haploid size of pop1
3) haploid size of pop2 

4) tab/space delimited list of site to be used -> chromosome position

5) mygroups -> number of groups for the deleted jackknife 

"""
def RAB(couplesoffreqs):
	# mycouplesoffreq=[(frq_sitei_popA, frq_sitei_popB), (frq_sitei+1_popA, frq_sitei+1_popB), (frq_sitei+2_popA, frq_sitei+2_popB), () ...]
	LAnotB=0; LBnotA=0
	for site in couplesoffreqs:
		freqA=site[0]; freqB=site[1] 
		LAnotB+=freqA*(1-freqB)
		LBnotA+=freqB*(1-freqA)
	RAB=LAnotB/float(LBnotA)
	return RAB 

def R2AB(couplesofderivedcounts, haploidsizepop1, haploidsizepop2):
	# mycouples=[ (counts_sitei_popA, counts_sitei_popB ), (counts_sitei+1_popA, counts_sitei+1_popB ), (counts_sitei+2_popA, counts_sitei+2_popB).... ]
	L2AnotB=0; L2BnotA=0	
	niA=haploidsizepop1; niB=haploidsizepop2 
	for site in couplesofderivedcounts: 
		diA=site[0];  diB=site[1]
		L2AnotB+= ( (diA*(1-diA) ) / (niA*(niA-1) )  ) * (1- (diB*(1-diB)) / (niB*(niB-1)) )
		L2BnotA+= ( (diB*(1-diB) ) / (niB*(niB-1) )  ) * (1- (diA*(1-diA)) / (niA*(niA-1)) )
	R2AB= L2AnotB/float(L2BnotA)
	return R2AB 

def deleted_mj_jackknife(theta,jackkpairs ) :
	"""REF http://www.soziologie.uni-halle.de/langer/buecher/mehrebenen/literatur/delete-m.pdf  
	theta is the estimation with all the observations 
	jackkpairs is a list of tuples (nb_of_deleted_observations, estimation_after_deletions_of_nb_obs)  """ 
        sum_a=0; sum_b=0
        mjs=[int(i[0]) for i in jackkpairs]
 	if len(mjs) >1: 
		n=sum(mjs) #TO ADD check if n equal total as in input   or send error message
        	g=len(mjs)
	        for i in jackkpairs:
        	        mj=i[0]; thetajstar=i[1]
                	a=(1-mj/float(n))*thetajstar; sum_a+=a
	        #print sum_a
        	thetajmj=g*theta-sum_a
	
	        for i in jackkpairs:
        		mj=i[0]; thetajstar=i[1]; hj=n/mj;
                	pseudothetajstar= hj*theta-(hj-1)*thetajstar;
 	               	b=(1/float(hj-1))*(pseudothetajstar-thetajmj)**2; sum_b+=b
        	#print sum_b
	        sigmasquarejmj=1/float(g)  *sum_b 
	        return thetajmj, sigmasquarejmj
	else: return "only one group" 

def constrained_sum_sample_pos(n, total):
    dividers = sorted(random.sample(xrange(1, total), n - 1))
    return [a - b for a, b in zip(dividers + [total], [0] + dividers)]

# /////////////////////////////////////////////

# ***** list of loci to be used *************
mysitestobeused={}

for line in open(mylistofsites, 'r').readlines():
	x=line.rstrip().split()
	chromosome=x[0]; position=x[1]
	if not chromosome in mysitestobeused: 
		mysitestobeused[chromosome]=[]
	mysitestobeused[chromosome].append(position)
	
# ***** derived allele counts  *************
mycouplesoffreqs=[]
mycouplesofcounts=[]

countloci=0
for line in open(myinputderivedallelecounts, 'r').readlines():
	myvec=line.rstrip().split()
	chrom=myvec[0]; pos=myvec[1]
	if pos in mysitestobeused[chrom]: 
		countloci+=1 
		f= [float(x) for x in line.rstrip().split()] 
       		mycouplesoffreqs.append( ( f[2]/float(mypop1haploidsize), f[3]/float(mypop2haploidsize) ) )
        	mycouplesofcounts.append(( f[2], f[3]) )

#print mycouplesofcounts
#print mycouplesoffreqs

# ***** RAB, R2AB, Jackknife   *************

n=len(mycouplesoffreqs) # number of observations 
g=mygroups     # number of groups in which you want to subdivide the  observations 
mjs=constrained_sum_sample_pos(g, n ) # subdivide the n observations in  g groups of sizes j in mjs 

jackkpairsRAB=[]
jackkpairsR2AB=[]

start=0 
for j in mjs:
	#print "/////////////////", j , start 
	stop=start+j
	temp_freqs=list(mycouplesoffreqs); # make a copy of the allele frequencies list as  temporary list 
	temp_freqs[start:stop]=[]; # delete start:j observations 

	temp_counts=mycouplesofcounts # make a copy of the allele counts list as  temporary list
	temp_counts[start:stop]=[] # delete start:j observations


	if len(temp_freqs)> 2:

		rabj = RAB(temp_freqs) # estimate RAB from the  subset of values 
		jackkpairsRAB.append((j, rabj))   # store the RAB value 

		r2abj = R2AB(temp_counts, mypop1haploidsize, mypop2haploidsize) # estimate R2AB from the  subset of values 
                jackkpairsR2AB.append((j, r2abj))   # store the RAB value

	else: 
		pass 
	start=stop
	#print "/////////////////", stop 

#listofjackkpairsRAB=[x[1] for x in jackkpairsRAB]


myRAB= RAB(mycouplesoffreqs); 
myjackrab=deleted_mj_jackknife(myRAB ,jackkpairsRAB )

myR2AB= R2AB(mycouplesofcounts, mypop1haploidsize, mypop2haploidsize);
myjackr2ab=deleted_mj_jackknife(myR2AB ,jackkpairsR2AB )


# ***** print results   *************



results=[str(x) for x in [ "RAB\t%s" %(myRAB), "jRAB\t%s" %(myjackrab[0]) , "sigmasquare\t%s" % (myjackrab[1]),  "R2AB\t%s" %(myR2AB), "jR2AB\t%s" %(myjackrab[0]) , "sigmasquareR2AB\t%s" % (myjackr2ab[1])]]

print "\n".join(results) 

