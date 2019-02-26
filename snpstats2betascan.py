import sys 
import argparse
import re 



def writeline(countA, countB, diploidTotal, imputeInfo, imputeInfoThreshold):
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

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Name of input snpstat file ",type=str,required=True)
    parser.add_argument("-i", help="impute info threshold",type=float, required= True)
    parser.add_argument("-o", help="pathto output file",type=str,required=True)
    args = parser.parse_args()
    output = open(args.o,'w')

    for line in open(args.f,'r'): 
	dataInLine=snpstatLineParser(line)
	if dataInLine is not None: 
		argsforwriteline=dataInLine[1:]+[args.i]
		resultscounts=writeline(*argsforwriteline)
		if resultscounts is not None: 
			results=[dataInLine[0]] + resultscounts+["\n"]
			#print " ".join(map(str,results )) 
                        output.write(" ".join(map(str,results )))


if __name__ == "__main__":
    main() 


