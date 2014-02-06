#! /usr/local/bin/python

import sys
import tarfile
from numpy import zeros, linalg

#dimensions=[ 161, 193]
#boxes=[ 19, 10 ]
boxes=[(30,45,40),(60,45,40)]
dimensions=[225]
dimeguess=.25

def getbox(dime,approx) :
    prevdiff=1e100
    for i in range(100) :
        diff=abs(approx-(1+32*i)*dime)
        if diff > prevdiff :
            boxsize=(1+32*(i-1))*dime
            break
        else :
            prevdiff=diff
    return boxsize,1+32*(i-1)

def getdime(box,approxdime) :
    prevdiff=1e100
    for i in range(100) :
        diff=abs(approxdime-box/(1+32.*i))
        if diff > prevdiff :
            dime=1+32.*(i-1)
            break
        else : prevdiff=diff
    return int(dime)

def usage():
    print "usage: %s pqr-archive.tar topol.itp" % sys.argv[0]
    sys.exit()
def resid(l) :
    return int(l[21:26])

def CoM( coordfile ) :
    RalCoords=[0,0,0]
    RalsumMass=0
    RapCoords=[0,0,0]
    RapsumMass=0
    cdcoords=0
    necoords=0
    for line in coordfile.readlines() :
        if "ATOM" in line and "DUM" not in line :
            coords=line.split()
            atom=coords[2]
            if len(atom) == 4 :
                atom = atom[1:] + atom[0]
            resname=coords[3]
            mass=atomname2gmx[resname,atom]
            if resid(line) < 101 :
                RalsumMass+=mass
                RalCoords[0] += float(coords[5])*mass
                RalCoords[1] += float(coords[6])*mass
                RalCoords[2] += float(coords[7])*mass
            elif resid(line) >= 100 :
                RapsumMass+=mass
                RapCoords[0] += float(coords[5])*mass
                RapCoords[1] += float(coords[6])*mass
                RapCoords[2] += float(coords[7])*mass
            if "CNC" in line :
                if "CD" in line :
                    cdcoords=line.split()
                    cdcoords=tuple( [ float( cdcoords[i] ) for i in (5,6,7) ] )
                elif "NE" in line :
                    necoords=line.split()
                    necoords=tuple( [ float( necoords[i] ) for i in (5,6,7) ] )
    RalCoords=tuple( [float(RalCoords[i])/RalsumMass for i in range(3)] )
    if RapsumMass == 0 : Coords=RalCoords
    else : 
        RapCoords=tuple( [float(RapCoords[i])/RapsumMass for i in range(3)] )
        Coords=tuple( [RapCoords[i]/2+RalCoords[i]/2 for i in range(3)] )

    return Coords, cdcoords, necoords

def inBox( coord, com, boxdimex, boxdimey, boxdimez ) :
    minx=com[0]-boxdimex/2.
    maxx=com[0]+boxdimex/2.
    miny=com[1]-boxdimey/2.
    maxy=com[1]+boxdimey/2.
    minz=com[2]-boxdimez/2.
    maxz=com[2]+boxdimez/2.
    mincoord=tuple([ minx, miny, minz ])
    maxcoord=tuple([ maxx, maxy, maxz ])
    for i in range(3) :
        if coord[i] < mincoord[i] or coord[i] > maxcoord[i] :
            print "(%.3f, %.3f, %.3f) outside of box! Box ranges from"%coord,mincoord,"to",maxcoord,
            if i == 0 : print ": x coord",coord[0]-mincoord[0]
            elif i == 1 : print ": y coord",coord[1]-mincoord[1]
            elif i == 2 : print ": z coord",cood[2]-mincoord[2]
            return False
            break
    

try :
    # topol.top or topol_A.itp
    topfile=open(sys.argv[2])
except :
    print "Could not open topology file." 
    usage()

atomname2gmx = {}
for line in topfile.readlines() :
    try :
        resname=line.split()[3]
        atom=line.split()[4]
        mass=float(line.split()[7])

        atomname2gmx[ resname, atom ] = mass
    except IndexError : pass # means this is not a line for an atom
    except ValueError : pass # means column 7 is not a mass
topfile.close()

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

for member in archive :
    pqrfile = archive.extractfile( member )
    framename = pqrfile.name
    #print framename,
    CoMcoords, CDcoords, NEcoords=CoM(pqrfile)
    #print "CoM: %.3f %.3f %.3f"%(CoMcoords),
    pqrfile.close()
    
    CDdist=tuple( [CoMcoords[i]-CDcoords[i] for i in range(3) ] )
    NEdist=tuple( [CoMcoords[i]-NEcoords[i] for i in range(3) ] )
    #print "CD: %.3f, NE: %.3f" %(linalg.norm(CDdist),linalg.norm(NEdist))
    if linalg.norm( NEdist ) > linalg.norm( CDdist ) : startpoint=NEcoords
    else : startpoint=CDcoords
    
    # vector from furthest atom to CoM
    vector=tuple( [ CoMcoords[i]-startpoint[i] for i in range(3) ] )
    normvector=vector/linalg.norm(vector)

    # replace PQRNAME with framename
    newapbslines_base = apbslines.replace("PQRNAME",framename)

    # for each box size, do
    # for each dimevalue, replace DIME and OUTPUT
    # OUTPUT is a combination of reference input file, pqrname, and dimevalue
    apbsoutname = "%s" % framename
    apbsoutname = apbsoutname.replace( ".pqr", ".in" )
    apbsout = open( apbsoutname, "w" )

    bcname = "%s" % ( framename )
    bcname = bcname.replace(".pqr","")
    
    read = "read\n\tmol pqr %s\n\tpot dx %s.bc.dx\nend\n\n" \
        % ( framename, bcname )

    apbsout.write( read )
    
    for box in boxes :
            for dimevalue in dimensions :
                # Center is the CoM -or-, (box/2-1) angstroms along the vector towards the furthest of CD and NE
                # Going to start at the furthest atom away to make sure everything is in the box
                #print linalg.norm(vector), vector
                #if linalg.norm(vector) <= (box/2*.9) : center=CoMcoords
                #else : center=tuple( [ startpoint[i]+(box/2*.9)*normvector[i] for i in range(3) ] )
                #print "CD coords: %.3f %.3f %.3f,"%CDcoords, "NE coords: %.3f %.3f %.3f," %NEcoords,
                #print "Center: %.3f %.3f %.3f" %center

                center=CoMcoords
                inBox( CDcoords, center, box[0],box[1],box[2] )
                inBox( NEcoords, center, box[0],box[1],box[2] )

                
                # replace CDCOORDS with new center coords
                newapbslines = newapbslines_base.replace( "CDCOORDS", "%.3f %.3f %.3f" %center )

                txtoutname = "%s_%s_%d" % (bcname, box[0], dimevalue)


                dimestring = "%d %d %d" %(getdime(box[0],dimeguess),getdime(box[1],dimeguess),getdime(box[2],dimeguess))
                newapbslines = newapbslines.replace("DIME",dimestring)

                boxstring="%.3f %.3f %.3f" %box
                newapbslines_0 = newapbslines.replace("FOCUS",boxstring)

                
                newapbslines_0 = newapbslines_0.replace( "OUTPUT", "%s_1" %txtoutname )


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



