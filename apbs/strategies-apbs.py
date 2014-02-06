#!/usr/bin/python
# Generate APBS input files (from a reference), by modification of four things:
# 1. for each pqr file in the tar.bz2 file, replace PQRNAME with the pqr name
# 2. for each pqr file in the tar.bz2 file, replace CDCOORDS with the
#   coordinates of the delta carbon of CNC
# 3. replace DIME with 161 and 193
# 4. replace OUTPUT with "%s_%d.dx" % ( framename, dimevalue )

import sys
import tarfile
from copy import copy

#dimevalues = [ 32*c+1 for c in range(3,7) ] # 97, 129, 161, 193
dimevalues = [ 32*c+1 for c in range(5,7) ] # 161, 193
boxvalues = [ i for i in range(10,20,9) ] # 10, 19

def usage(): 
    print "usage: %s pqr-archive.tar" % sys.argv[0]
    sys.exit()

try:
    inputFileName = "/Users/ritchie/Utilities/apbs/strategies.in" # name of APBS .in file
    BCFileName = "/Users/ritchie/Utilities/apbs/BC.in"
    archiveName = sys.argv[1]   # name of pqr archive

    # read the archive
    archive = tarfile.open( archiveName )

    # read the APBS input file
    apbsinfile = open( inputFileName )
    apbslines = apbsinfile.read()
    apbsinfile.close()

    # read the BC input file
    bcinfile = open( BCFileName )
    bclines = bcinfile.read()
    bcinfile.close()

except IndexError, TypeError :
    usage()
    sys.exit()
except IOError :
    print "could not open %s or %s" % ( inputFileName, archiveName )
    usage()
except :
    print "some error"
    usage()

# Loop over archive members
for member in archive:
    pqrfile = archive.extractfile( member )
    framename = pqrfile.name
    
    # get CNC-CD coordinates
    for line in pqrfile.readlines() :
        if "CNC" in line :
            if "CD" in line :
                cdcoords = line.split()
                cdcoords = tuple( [ float( cdcoords[i] ) for i in (5,6,7) ] )
                break
    pqrfile.close()

    # replace CDCOORDS with the newly read coordinates
    newapbslines_base = apbslines.replace( "CDCOORDS", "% .3f % .3f % .3f"  % cdcoords )
    
    # replace PQRNAME with framename
    newapbslines_base = newapbslines_base.replace( "PQRNAME", framename )                                      

    # for each box size, do
    # for each dimevalue, replace DIME and OUTPUT
    # OUTPUT is a combination of reference input file, pqrname, and dimevalue
    apbsoutname = "%s" % framename
    apbsoutname = apbsoutname.replace( ".pqr", ".in" )
    apbsout = open( apbsoutname, "w" )

    bcname = "%s" % ( framename )
    bcname = bcname.replace( ".pqr", "" )

    read = "read\n\tmol pqr %s\n\tpot dx %s.bc.dx\nend\n\n" \
        % ( framename, bcname )

    apbsout.write( read )

    for box in boxvalues :
        for dimevalue in dimevalues :
            txtoutname = "%s_%s_%d" % ( bcname, box, dimevalue )

            newapbslines = newapbslines_base.replace( "OUTPUT", txtoutname )

            dimestring = "%d %d %d" % ( dimevalue, dimevalue, dimevalue )
            newapbslines = newapbslines.replace( "DIME", dimestring )

            boxstring = "%d %d %d" % ( box, box, box )
            newapbslines_0 = newapbslines.replace( "FOCUS", boxstring )

            apbsout.write( newapbslines_0 )
            
    apbsout.write( "quit" )
    apbsout.close()

    # Build boundary condition input file
    bcoutname = "%s" % framename
    bcoutname = bcoutname.replace( ".pqr", ".bc.in" )
    bcout = open( bcoutname, "w" )

    # replace PQRNAME with framename
    newbclines = bclines.replace( "PQRNAME", framename )

    # OUTPUT is the same as before, but will print with a .bc suffix
    newbclines_base = newbclines.replace( "OUTPUT", bcname )

    bcout.write( newbclines_base )
    bcout.close()

archive.close()
