import sys 
import argparse
import re 
import numpy as np



def MinorAlleleCounts(countA, countB, diploidTotal, imputeInfo, imputeInfoThreshold):
    """ 1. establish which is the minor allele and approximate to the haploid total
	2. returns results if impute_info>= impute_info_threshold """	 
    haploidTotal=int(diploidTotal)*2 
    majorAllele=max(float(countA), float(countB) )
    minAllele=int(haploidTotal-majorAllele)
    if not imputeInfo=='NA': #apprently in monomorphic sites 
        if  float(imputeInfo)>=float(imputeInfoThreshold): 
            return [minAllele, haploidTotal]

def snpstatLineParser(line): 
    if not re.search( '#|position', line ):
        x=line.rstrip().split(" ") 
        position=x[3]; imputeInfo=x[17] 
	countA=x[9];  countB=x[10]; diploidTotal=x[25]
        return [position, countA, countB,  diploidTotal, imputeInfo]  
    	#return(position)    


def snpsPerWindow(listofpositions, coreSnpPos, window):
#"""modified from https://stackoverflow.com/questions/54237254/count-values-in-overlapping-sliding-windows-in-python"""
    bin_starts = coreSnpPos-window/float(2) #np.arange(start, end+1-window, step)
    bin_ends = bin_starts + window
    last_index = np.searchsorted(listofpositions, bin_ends, side='right')
    first_index = np.searchsorted(listofpositions, bin_starts, side='left')
    return  [int(bin_starts), int(bin_ends), last_index - first_index]  

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Name of input snpstat file ",type=str,required=True)
    parser.add_argument("-i", help="impute info threshold",type=float, required= True)
    parser.add_argument("-o", help="pathto output file",type=str,required=True)
    parser.add_argument("-maf", help="Minor Allele Frequency of the core SNPs",type=float,required=True)
    parser.add_argument("-w", help="window size around core SNPs",type=int,required=True)
    args = parser.parse_args()

    output = open(args.o,'w')

    allPositions=[]
    corePositions=[]

    for line in open(args.f,'r'): 
	dataInLine=snpstatLineParser(line)
	if dataInLine is not None: 
		allPositions.append(int(dataInLine[0])) 
		argsforMinorAlleleCounts=dataInLine[1:]+[args.i]
		resultscounts=MinorAlleleCounts(*argsforMinorAlleleCounts)
		if resultscounts is not None: 
			#print (resultscounts[1] , float(resultscounts[0]) 
	                if resultscounts[0]/float(resultscounts[1]) >= args.maf: 
                            #print (( resultscounts[1], resultscounts[0]) ) 
			    corePositions.append(int(dataInLine[0]))
		#	results=[dataInLine[0]] + resultscounts+["\n"]
			#print " ".join(map(str,results )) 
                 #       output.write(" ".join(map(str,results )))
    
    output.write( " ".join (["impute_info", "window_bp", "maf", "positionStart", "positionEnd",  "nSNPs", "\n"])) 
    if len(corePositions)>0: 
	    for i in corePositions:
		results= [args.i, args.w, args.maf] + snpsPerWindow(allPositions,int(i) , args.w )  +['\n']
		output.write( " ".join(map(str,results )))	
    else: print ('not enough  SNPs with maf > %s and impute info > %s ' % (args.i,  args.maf) )		
    

if __name__ == "__main__":
    main() 


