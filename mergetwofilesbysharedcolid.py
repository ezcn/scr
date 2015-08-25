import pandas as pd
import sys 

mysharedcolumnid=sys.argv[1]
myfirstfile=sys.argv[2]
mysecondfile=sys.argv[3]

file_one = pd.read_csv(myfirstfile,sep='[ \t]', keep_default_na=False, na_values=[""], engine='python')
file_two = pd.read_csv(mysecondfile, sep='[ \t]', keep_default_na=False, na_values=[""], engine='python')


merged_inner = pd.merge(left=file_one,right=file_two, on=mysharedcolumnid)
print merged_inner 

