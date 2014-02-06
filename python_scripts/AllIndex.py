#! /usr/bin/python

import sys

try:
    file=open(sys.argv[1])
except:
    print "select input file"
    sys.exit()

pdb=file.readlines()[:-1]
file.close()
def inCoord(line) :
        inline=False
        try :
            float(line.split()[-1])
            inline=True
        except :
            return False
        try :
            float(line.split()[-2])
            inline=True
        except :
            return False
        try :
            float(line.split()[-3])
            inline=True
        except :
            return False
        try :
            float(line.split()[-4])
            inline=True
        except :
            return False
        if inline :
            return True
print "[ complex ]"
for line in pdb:
    if inCoord(line) :
        ndx=line[15:20]
        print ndx

