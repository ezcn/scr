import sys 
import argparse
import re 



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

def snpsPerWindow (snpPositions, coreSnp, windowsize): 
    positionStart=coreSnp-windowsize/float(2) #min(snpPositions)
    positionEnd=coreSnp+windowsize/float(2)#max(snpPositions)
    return [int(positionStart), int(positionEnd), len([x for x in snpPositions if positionStart <= x <= positionEnd])] 

   

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
		results= [args.i, args.w, args.maf]+snpsPerWindow(allPositions, int(i), args.w) +['\n']
		output.write( " ".join(map(str,results )))	
    else: print ('not enough  SNPs with maf > %s and impute info > %s ' % (args.i,  args.maf) )		
    

if __name__ == "__main__":
    main() 


