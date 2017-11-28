import sys, re, numpy, pandas , random 


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

#''''''''''''''''''''''''''''''''''''''''''''''

data_set=[random.uniform(20, 100) for i in range( 100) ]  
#print data_set 

print numpy.mean(data_set) 

del_obs=constrained_sum_sample_pos(4, 100)
print del_obs

start=0
jackkdata=[]
for j in del_obs:
	print "/////////////////", j , start 
	stop=start+j
	temp_data=list(data_set); # make a copy of the data  list as  temporary list 
	temp_data[start:stop]=[]; # delete start:j observations 
	print 'len temp data' , len(temp_data) 

	if len(temp_data)> 2:
		mystat=numpy.mean(temp_data) # estimate the stats  from the  subset of values 
		jackkdata.append((j, mystat)) # store the number of deleted observations and the  stats  

	else: 
		pass 
	start=stop

print jackkdata
print deleted_mj_jackknife(numpy.mean(data_set) ,jackkdata ) 
