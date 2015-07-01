#! /usr/bin/env python
# writedummies.py
# put dummy atoms between CD and NE of the nitrile
# dummy atoms include CD and NE coordinates

import sys
from numpy import array, dot, sqrt
    
usage = "python[2] read.py pdb.pqr\n";
nsteps = 10

try :
    pqrfile = sys.argv[1]
except IndexError :
    print "*** Syntax error: got %d arguments, expected 2." % len(sys.argv)
    print "%s\n" % usage
    sys.exit(2)

# now find CNC NE and CNC CD coordinates in the pdb file
file = open( pqrfile )
pqrlines = file.readlines()
file.close()
yes1 = yes2 = False
for line in pqrlines :
    if "CNC" in line :
        if "CD" in line :
            cdcoords = float( line[30:38] ), float( line[38:46] ), float( line[46:54] ) 
            cdcoords = array(cdcoords)
            yes1 = True
        if "NE" in line :
            necoords = float( line[30:38] ), float( line[38:46] ), float( line[46:54] ) 
            necoords = array(necoords)
            yes2 = True

    # if we have found CD and NE, then write nsteps dummy atoms
    #ATOM     11  H   SER     2      41.040  36.180  22.490  0.3454   0.600
    #01234567890123456789012345678901234567890123456789012345678901234567890 
    #0         1         2         3         4         5         6         7
    if yes1*yes2 : 
        dummylines = ""
        stepvector = (necoords - cdcoords )/nsteps
        for i in range(0,nsteps+1):
            dummycoords = cdcoords + stepvector*i
            dummyline = "ATOM      0 DUM  CNC     0    % 8.5f% 8.5f% 8.5f  0.0000   0.000\n" % tuple(dummycoords)
            dummylines += dummyline
        yes1 = yes2 = False        

# write the new pqr with the dummy atoms at the beginning
print dummylines.strip()
for line in pqrlines : print line.strip()
