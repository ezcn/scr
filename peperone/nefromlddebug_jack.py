import sys, re, numpy, pandas 

"""
- McEvoy BP, Powell JE, Goddard ME, Visscher PM (2011) Genome research 21: 821-829 
- Hayes BJ, Visscher PM, McPartlan HC, Goddard ME (2003) Genome Research 13: 635-643 

''''''''''''''''''''''''''''''''''''''''
input file: first line is ignored;  as many column as wanted ony first two taken into account: 
- first column:  distance in cM 
- second column: r2  

usage:   python Ne_LD.py infile n_haplo """

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

infile = sys.argv[1]
n_haplo = sys.argv[2] # number of chromosomes used to calculate ld 
min_nb_observ_in_bin = 2 # minimum nuber of observations in a bin 
maxdist_cM = 50  # 0.5 max distance in cM between marker pairs 
step_cM = 0.001   #0.005 increment in cM for bins 
minr2 = 0.0000001 # minimun r2 to be considered  
maxr2 = 0.9999 # maximum r2 to be considered


#''''''''''' functions 

def ne_estim (rsq, morg_dist): 
	"""McEvoy(2011)GenRes 
	rsq = average r^2 for marker pairs with similar genetic distance
	morg_dist = genetic distance in Morgan (average, similar) between marker pairs  
	Ne ~ 1/ (4c) * [(1/r2)-2]"""

	ne = 1/(4*morg_dist) * ( (1/float(rsq) ) -2 )
	return ne


def ne_estim_corr (rsq, morg_dist, n ): 
 	"""McEvoy(2011)GenRes 
        rsq = average r^2 for marker pairs with similar genetic distance
        morg_dist = genetic distance in Morgan (average, similar) between marker pairs  
        n = haploid sample size  (number of chromosomes) 
        Ne ~  1/ (4c) * [(1/(r2-(1/n)))-2]"""

	a = 1/(4*morg_dist)
        b = rsq - (1/float(n) )
        ne_c = a * ( (1/float(b) ) -2)
        return ne_c

def t_estim (morg_dist): 
	"""Hayes 2003 GenRes
	morg_dist = genetic distance in Morgan"""

	t = 1/(2*morg_dist)
	return t 

def constrained_sum_sample_pos(n, total):
	dividers = sorted(random.sample(xrange(1, total), n - 1))
	return [a - b for a, b in zip(dividers + [total], [0] + dividers)]	

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
        		mj=i[0]; thetajstar=i[1]; hj=n/float(mj);
                	pseudothetajstar= hj*theta-(hj-1)*thetajstar;
 	               	b=(1/float(hj-1))*(pseudothetajstar-thetajmj)**2; sum_b+=b
        	#print sum_b
	        sigmasquarejmj=1/float(g)  *sum_b 
	        return thetajmj, sigmasquarejmj
	else: return "only one group" 


#''''''''''''  define bins of distance in cM 

bins = [i for i in numpy.arange (0, maxdist_cM, step_cM )]
binsnames = [i for i in range(1, len(bins) ) ]
#binsnames =map(chr, range(65, 65+len(bins)-1 ))   

#''''''''''   bin the input file 
raw_data = {'distance':[], 'r2':[]}

with open (infile, "r") as f: 
	next(f) 
	for line in f: 	
		x = line.split()
		indist = float(x[0]); inr2=float(x[1])
		if minr2 <= inr2 <= maxr2:  
			raw_data['distance'].append(indist) 
			raw_data['r2'].append(inr2)

df = pandas.DataFrame(raw_data, columns = ['distance', 'r2'])
df['categories'] = pandas.cut(df['distance'], bins,  labels=binsnames, include_lowest=True)

cat_counts=pandas.value_counts(df['categories'])

#with pandas.option_context('display.max_rows', None, 'display.max_columns', 5):
#	print(df) 	
 

#'''''''''''  iterate over bins of the df and print results 

title=['bin','nb_marker_pairs' ,  'average_r2', 'average_distance_cM' , 'time_gen', 'ne' , 'ne_corr']
print "\t".join(title)

for b in binsnames:
	r2set = df.loc[df['categories']==b]['r2'] # select all r2 values in the bin 
	distset = df.loc[df['categories']==b]['distance']

	if len(r2set) >= min_nb_observ_in_bin: 
		meanr2 = numpy.mean(r2set) # average r2 in bin 
		bin_distance = numpy.mean(distset)/100.0 ## average distance in Morgan in the bin		

		if bin_distance >0 :  # if genetic distance is 0M estimates are not possibles  
			time = t_estim( bin_distance) 
			ne_corr = ne_estim_corr(meanr2, bin_distance , n_haplo)  
			ne = ne_estim(meanr2, bin_distance)  

			res=[b, len(r2set) , meanr2, bin_distance, time, ne, ne_corr ]
        		print "\t".join(map(str, res))
