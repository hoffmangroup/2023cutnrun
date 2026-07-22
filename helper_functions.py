import os
import subprocess

# finds the number of reads of a bam file. it would be an improvement if 
# we could figure out a way to count reads of a bam file purely using python, 
# and not use a separate bash command to do this. 
def BamReads(fname):    
    p = subprocess.run("sambamba flagstat " + fname, shell=True, stdout=subprocess.PIPE)
    mapped_counts = int(p.stdout.decode("utf-8").split("\n")[4].split()[0])
    printsep(fname + " reads: " + str(mapped_counts))
    return mapped_counts

# just an easier way to separate print statements. 
def printsep(input = None):
    print('')
    if(input == None):
        print('-------------------------------------------------------------------------')
    else:
        print(input)
    print('')
    