#! /usr/local/bin/python

import os
import sys
from numpy import array, linalg, polyfit, poly1d, polyder, min, arange, sum, product
from math import pi
from optparse import OptionParser
from multiprocessing import Pool

parser=OptionParser()
parser.add_option("-f",\
                  "--directory",\
                  dest="directory",\
                  help="directory name"\
                  )
parser.add_option("-t",\
                  "--top",\
                  dest="topfile",\
                  help=".top or .itp file (for mass information)"\
                  )
parser.add_option("-d",\
                  "--docked",\
                  dest="dtemp",\
                  default="/Users/ritchie/Desktop/tmp_apbs/make_inputs/template.in",\
                  help="Template.in for docked complex"\
                  )
parser.add_option("-m",\
                  "--monomer",\
                  dest="mtemp",\
                  default="/Users/ritchie/Desktop/tmp_apbs/make_inputs/monomer_template.in",\
                  help="Template.in for monomer"\
                  )

option,args=parser.parse_args()

dtemp=option.dtemp
mtemp=option.mtemp
directory=option.directory
topfile=option.topfile

directory="/Volumes/Jotunheim/ritchie/MD/2D_Umbrella/0-E30D_K31E+N27C/"
topfile="/Volumes/Jotunheim/ritchie/MD/2D_Umbrella/0-E30D_K31E+N27C/solv_ion-0-0_A.itp"

def getBondVector( structurefilelines ) :
    try :
        for line in structurefilelines:
            if "CNC" in line and "CD" in line :
                cdcoords=line.split()
                cdcoords=array([float(cdcoords[i]) for i in (5,6,7)])
            if "CNC" in line and "NE" in line :
                necoords=line.split()
                necoords=array([float(necoords[i]) for i in (5,6,7)])
        BondVector = necoords-cdcoords
    except :
        print "ERROR: Empty pqr file!"
        sys.exit()
    return BondVector, cdcoords, necoords

def resid(l) :  return int(l.split()[4])

def getdime(box,approxdime) :
    dime=[]
    for n in range(len(box)) :
        prevdiff=1e100
        for i in range(100) :
            diff=abs(approxdime-box[n]/(1+32.*i))
            if diff > prevdiff :
                dime.append(1+32*(i-1))
                break
            else : prevdiff=diff
    return dime

def box_coords( coord, dimensions ) :
    return coord-dimensions/2, coord+dimensions/2
    

def inBox( coord, com, boxdime ) :
    mincoord,maxcoord = box_coords( com, boxdime )
    for i in range(3) :
        if coord[i] < mincoord[i] or coord[i] > maxcoord[i] :
            print "(%.3f, %.3f, %.3f) outside of box! Box ranges from"%(coord[0],coord[1],coord[2]),mincoord,"to",maxcoord
            inbox = False
            sys.exit()
        else : inbox = True
    return inbox

def CoM( coordlines, min_resid=0, max_resid=1e100 ) :
    Coords=array([0.,0.,0.])
    sumMass=0
    for line in coordlines:
        if "ATOM" in line and "DUM" not in line and min_resid <= resid(line) <= max_resid :
            coords=line.split()
            atom=coords[2]

            resname=coords[3]
            mass=atomname2gmx[resname,atom]
            sumMass+=mass
            Coords[0] += float(coords[5])*mass
            Coords[1] += float(coords[6])*mass
            Coords[2] += float(coords[7])*mass
            if "CNC" in line :
                if "CD" in line :
                    cdcoords=line.split()
                    cdcoords=array( [ float( cdcoords[i] ) for i in (5,6,7) ] )
                elif "NE" in line :
                    necoords=line.split()
                    necoords=array( [ float( necoords[i] ) for i in (5,6,7) ] )
    Coords=array( [float(Coords[i])/sumMass for i in range(3)] )
    try :
        return Coords, cdcoords, necoords, necoords-cdcoords
    except :
        return Coords

try : top=open(topfile)
except :
    print "Could not open topology file, %s" %topfile
    sys.exit()

def get_box_center( com, cdcoords, necoords, box_array ) :
    CDdist = linalg.norm( com - cdcoords )
    NEdist = linalg.norm( com - necoords )
    if NEdist > CDdist : start_point = necoords
    else : start_point = cdcoords
    
    translation_vector = com - start_point
    vector_mag = linalg.norm(translation_vector)
    normal_vector = translation_vector / vector_mag
    smallest_dimension = min(box_array)
    
    if vector_mag <= smallest_dimension/2*.9 : center = com
    else : center = start_point + normal_vector*box_array/2*.9

    if inBox( cdcoords, center, box_array ) and inBox( necoords, center, box_array ) :
        return center

def coordComp(a,b):
	same=True
	for i in xrange(len(a)) :
		if a[i] != b[i] :
			return False
		else :
			same=True
	return True


def make_inputs( (r1,r2) ) :
    print "Working on %s"%(os.path.join(directory,"%i-%i/"%(r1,r2)))
    for n in range(81) :
        pqrfile=os.path.join(directory,"apbs/%i-%i/%s_%i-%i-frame%i.pqr"%(r1,r2,name,r1,r2,n))
        pqr=open(pqrfile)
        pqrlines=pqr.readlines()
        pqr.close()
        SCoM,CD,NE,BV=CoM(pqrlines)
        SCoM10 = get_box_center( SCoM, CD, NE, b10 )
        SCoM19 = get_box_center( SCoM, CD, NE, b19 )
        RalCoM,CD,NE,BV=CoM(pqrlines,max_resid=100)
        RCoM10 = get_box_center( RalCoM, CD, NE, b10 )
        RCoM19 = get_box_center( RalCoM, CD, NE, b19 )

        if not coordComp(RalCoM,SCoM) :
            RapCoM=CoM(pqrlines,min_resid=101)
            interface = (RapCoM+RalCoM)/2
            templatelines = dockedlines
        else :
            interface = SCoM
            template = monomerlines
        
        outname=pqrfile.replace(".pqr","")
        inname=pqrfile.replace(".pqr","_2.in")
        inname2=inname.replace(".in","_mg10.in")
        inname3=inname.replace(".in","_mg19.in")

        outfile=open(inname,"w")
        outfile2=open(inname2,"w")
        outfile3=open(inname3,"w")
        first=True
        second=False
        third=False
        for line in templatelines :
            if "mg-auto" in line and "#" in line and first and " 10 " in line :
                first=False
                second=True
                third=False
                outfile2.write("read\n\tmol pqr %s\nend\n\n"%pqrfile)
                outfile.write("quit\n")
            elif "mg-auto" in line and "#" in line and second and " 19 " in line :
                first=False
                second=False
                third=True
                outfile3.write("read\n\tmol pqr %s\nend\n\n"%pqrfile)
                outfile2.write("quit\n")
            if first :
                outfile.write( line.replace("OUTPUT",outname).\
                   replace("PQRNAME",pqrfile).\
                   replace("SCOM10","%.4f %.4f %.4f" %tuple(SCoM10)).\
                   replace("SCOM19","%.4f %.4f %.4f" %tuple(SCoM19)).\
                   replace("SCOM","%.4f %.4f %.4f" %tuple(interface)).\
                   replace("RCOM10","%.4f %.4f %.4f" %tuple(RCoM10)).\
                   replace("RCOM19","%.4f %.4f %.4f" %tuple(RCoM19)).\
                   replace("CDCOORDS","%.4f %.4f %.4f" %tuple(CD)))
            elif second :
                outfile2.write( line.replace("OUTPUT",outname).\
                   replace("PQRNAME",pqrfile).\
                   replace("SCOM10","%.4f %.4f %.4f" %tuple(SCoM10)).\
                   replace("SCOM19","%.4f %.4f %.4f" %tuple(SCoM19)).\
                   replace("SCOM","%.4f %.4f %.4f" %tuple(interface)).\
                   replace("RCOM10","%.4f %.4f %.4f" %tuple(RCoM10)).\
                   replace("RCOM19","%.4f %.4f %.4f" %tuple(RCoM19)).\
                   replace("CDCOORDS","%.4f %.4f %.4f" %tuple(CD)))
            elif third :
                outfile3.write( line.replace("OUTPUT",outname).\
                   replace("PQRNAME",pqrfile).\
                   replace("SCOM10","%.4f %.4f %.4f" %tuple(SCoM10)).\
                   replace("SCOM19","%.4f %.4f %.4f" %tuple(SCoM19)).\
                   replace("SCOM","%.4f %.4f %.4f" %tuple(interface)).\
                   replace("RCOM10","%.4f %.4f %.4f" %tuple(RCoM10)).\
                   replace("RCOM19","%.4f %.4f %.4f" %tuple(RCoM19)).\
                   replace("CDCOORDS","%.4f %.4f %.4f" %tuple(CD)))
        outfile.close()
        outfile2.close()
        outfile3.close()


dockedt=open(dtemp)
dockedlines=dockedt.readlines()
dockedt.close()

monomert=open(mtemp)
monomerlines=monomert.readlines()
monomert.close()

directory=os.path.abspath(directory)
path=directory.split('/')
name=path[-1]
b19=array([19.,19.,19.])
b10=array([10.,10.,10.])
d193=array([193.,193.,193.])
b30=array([30.,45.,40.])
b60=array([60.,45.,40.])
atomname2gmx = {}
for line in top.readlines() :
    try :
        resname=line.split()[3]
        atom=line.split()[4]
        mass=float(line.split()[7])

        atomname2gmx[ resname, atom ] = mass
    except IndexError : pass # means this is not a line for an atom
    except ValueError : pass # means column 7 is not a mass
top.close()

if __name__ == '__main__':
    pool = Pool(8)
    TASKS = [(i,j) for i in range(0,359,30) for j in range(0,359,30)]
    pool.map(func=make_inputs,iterable=TASKS)
