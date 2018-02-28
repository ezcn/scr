import sys, re, numpy, pandas 

"""
- McEvoy BP, Powell JE, Goddard ME, Visscher PM (2011) Genome research 21: 821-829 
- Hayes BJ, Visscher PM, McPartlan HC, Goddard ME (2003) Genome Research 13: 635-643 

''''''''''''''''''''''''''''''''''''''''
input file: first line is ignored;  as many column as wanted ony first two taken into account: 
- first column:  distance in cM 
- second column: r2  

usage:   python Ne_LD.py infile n_haplo keyword """

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

infile = sys.argv[1]
n_haplo = sys.argv[2] # number of chromosomes used to calculate ld 
keyword = sys.argv[3:] # something to identify the samples 
min_nb_observ_in_bin = 20 # minimum nuber of observations in a bin 
maxdist_cM =2        # 0.25 for 1kg   # 0.5 max distance in cM between marker pairs 
step_cM = 0.01            # 0.01 for 1kg 	   #0.005 increment in cM for bins 
minr2 = 0.0001 # minimun r2 to be considered  
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
	

#''''''''''''  define bins of distance in cM 
maxdist_M = maxdist_cM/100.0 ; step_M = step_cM/100.0 # convert intervals in Morgan
bins = [i for i in numpy.arange (0, maxdist_M, step_M )]
binsnames = [i for i in range(1, len(bins) ) ]
#binsnames =map(chr, range(65, 65+len(bins)-1 ))   

#''''''''''   bin the input file 
raw_data = {'distance':[], 'r2':[]}

with open (infile, "r") as f: 
	next(f) 
	for line in f: 	
		x = line.split()
		indist = float(x[0])/100.0 ; inr2=float(x[1]) # store distance in Morgans 
		if minr2 <= inr2 <= maxr2:  
			raw_data['distance'].append(indist) 
			raw_data['r2'].append(inr2)

df = pandas.DataFrame(raw_data, columns = ['distance', 'r2'])
df['categories'] = pandas.cut(df['distance'], bins,  labels=binsnames, include_lowest=True)

cat_counts=pandas.value_counts(df['categories'])

#with pandas.option_context('display.max_rows', None, 'display.max_columns', 5):
#	print(df) 	
 

#'''''''''''  iterate over bins of the df and print results 

title=['bin','nb_marker_pairs' ,  'average_r2', 'average_distance_cM' , 'time_gen', 'ne' , 'ne_corr', 'species', 'clade', 'type']
print "\t".join(title)

for b in binsnames:
	r2set = df.loc[df['categories']==b]['r2'] # select all r2 values in the bin 
	distset = df.loc[df['categories']==b]['distance']
	#print "########" 
	#print distset 

	if len(r2set) >= min_nb_observ_in_bin: 
		meanr2 = numpy.mean(r2set) # average r2 in bin 
		#print '~~~~~~~~~~'
		#print numpy.mean(distset)  , ' mean dist' 
		bin_distance = numpy.mean(distset)   # average distance in Morgan in the bin		
		#print bin_distance,  ' bin dist'
		if bin_distance >0 :  # if genetic distance is 0M estimates are not possibles  
			time = t_estim( bin_distance) 
			ne_corr = ne_estim_corr(meanr2, bin_distance , n_haplo)  
			ne = ne_estim(meanr2, bin_distance)  

			res=[b, len(r2set) , meanr2, bin_distance*100, time, ne, ne_corr] + keyword
        		print "\t".join(map(str, res))
