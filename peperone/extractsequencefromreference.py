from pyfasta import Fasta 
import sys 

myfastaseq=sys.argv[1] 
mycoordinates=sys.argv[2:] 

header=[">"]+mycoordinates
print "_".join(map(str, header))
 
mychr=mycoordinates[0]; mystart=int(mycoordinates[1]) ; myend= int(mycoordinates[2]) 

f=Fasta(myfastaseq) 
print f.sequence({'chr': mychr , 'start': mystart, 'stop': myend}) 
