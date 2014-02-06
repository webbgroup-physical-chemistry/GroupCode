import sys
from numpy import array, linalg, polyfit, poly1d, polyder, min, arange, sum, product, ceil, floor
from math import pi


bctemplate="read \n\
    mol pqr PQRNAME \n\
end \n\
 \n\
# stage 4 \n\
# based on suggestions from Gernot Kierseritzky and stage4 described \n\
# in stages-logic.txt; changed to lpbe for boundary and convergence \n\
# issues when we go to finer grids, and spl2 because of fewer off- \n\
# the-grid warnings. \n\
elec \n\
    mg-manual \n\
 \n\
    dime 97 97 97 \n\
    glen 240.0 240.0 240.0 # 2.5 Angstroms \n\
 \n\
    gcent mol 1 \n\
    mol 1 \n\
    lpbe \n\
    bcfl sdh # faster than mdh and probably good enough if the coarse grid is large enough \n\
 \n\
    ion charge 1 conc 0.150 radius 2.0 \n\
    ion charge -1 conc 0.150 radius 2.0 \n\
    pdie 2.0 \n\
    sdie 78.0 \n\
    chgm spl2 \n\
    srfm mol \n\
    srad 1.4 \n\
    sdens 10.0 \n\
    temp 300.0 \n\
    calcenergy total \n\
    calcforce no \n\
    write pot dx OUTPUT.bc \n\
end \n\
 \n\
quit \n"

template="# stage2 \n\
# DIME \n\
# FOCUS \n\
elec \n\
    mg-manual \n\
 \n\
    glen FOCUS \n\
    dime DIME \n\
    gcent CDCOORDS \n\
    mol 1 \n\
    lpbe \n\
    bcfl map \n\
    usemap pot 1 \n\
 \n\
    ion charge 1 conc 0.150 radius 2.0 \n\
    ion charge -1 conc 0.150 radius 2.0 \n\
    pdie 2.0 \n\
    sdie 78.0 \n\
    chgm spl2 \n\
    srfm mol \n\
    srad 1.4 \n\
    sdens 10.0 \n\
    temp 300.0 \n\
    calcenergy total \n\
    calcforce no \n\
    write atompot flat OUTPUT \n\
end \n"

readdx="# read in stage1 calculation \n\
read \n\
    mol pqr PQRNAME \n\
    pot dx DXNAME \n\
end \n"


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

def make_bc_input( pqrname, dowrite=True ) :
    output=pqrname.replace(".pqr","")
    dxname="%s.bc.dx" %output
    newbc=bctemplate.replace("PQRNAME",pqrname).replace("OUTPUT",output)
    bcoutname=pqrname.replace(".pqr",".bc.in")
    if dowrite :
        bcout=open( bcoutname, "w" )
        bcout.write( newbc )
        bcout.close()
    return bcoutname,dxname

def make_apbs_input( pqrname, dxname, center_array, box_array, spacing_array, nlines, index=-1, first=False, dowrite=True ) :
    inputname=pqrname.replace(".pqr",".in")

    if index > -1 :
        output=pqrname.replace(".pqr","_%i_%i_%i" %(box_array[0],spacing_array[0],index))
    else :
        output=pqrname.replace(".pqr","_%i_%i" %(box_array[0],spacing_array[0]))
        
    if dowrite :
        if first :
            inputfile=open(inputname,"w")
            header=readdx.replace("PQRNAME",pqrname).replace("DXNAME",dxname)
            inputfile.write(header)
        else : inputfile=open(inputname,"a")

        try :
            completed=open('%s.txt'%output)
            completedlines=completed.readlines()[4:]
            nclines=len(completedlines)
            if nclines == nlines :
                print '%s.txt has been created, skipping...' %output
                completed.close()
                return "%s.txt" %(output), inputfile
            else : raise Exception()
        except :
            print "Setting up input file to include %s" %output
            boxstring="%.3f %.3f %.3f" %(box_array[0],box_array[1],box_array[2])
            dimestring="%d %d %d" %(spacing_array[0],spacing_array[1],spacing_array[2])
            coordstring="%.3f %.3f %.3f" %(center_array[0],center_array[1],center_array[2])
            newapbs=template.replace("DIME",dimestring).replace("FOCUS",boxstring).replace("OUTPUT",output).replace("CDCOORDS",coordstring)

            inputfile.write(newapbs)
    else : inputfile=""
        
    return "%s.txt" %(output), inputfile

def resid(l) :  return int(l[21:26])

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

def CoM( coordlines, atomname2gmx, min_resid=0, max_resid=1e100 ) :
    Coords=array([0.,0.,0.])
    sumMass=0
    for line in coordlines:
        if "ATOM" in line and "DUM" not in line and min_resid <= resid(line) <= max_resid :
            coords=line.split()
            atom=coords[2]
            if len(atom) == 4 :
                atom = atom[1:] + atom[0]
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
 
def parse_pot( pqrlines, potentialfile, bondlength, box_array, center, nlines, E30_resid, K31_resid ) :
    try :
        potfile=open(potentialfile)
        potlines=potfile.readlines()[4:]
        potfile.close()
    except :
        print "Error opening %s, exiting..." %potentialfile
        sys.exit()

    pot=[]
    pos=[]
    for i in arange(11) : # because only the first 11 atoms are our dummy atoms
        potential = float(potlines[i]) 
        if potential == 0 :
            error = True
        pot.append(potential)
        pos.append(bondlength*i/10.)
    pot,pos=array(pot),array(pos)
    
    if len(pos) == 0 or len(pot) == 0 :
        error = True

    # linear fit of middle 5 points
    linfield = -polyfit(pos[3:8],pot[3:8],1)[0]

    # 10th order polynomial fit of all 10 points at the midpoint
    polynomial = poly1d(polyfit(pos,pot,10))
    deriv = polyder(polynomial)
    polyfield = -deriv(0.5*bondlength)
    # 10th order polynomial fit at the inflection point
    seq = arange( 0.4*bondlength, 0.65*bondlength, 0.001 )
    inflfield = min(map(lambda x: -deriv(x), seq))

    charge,pcharge,ncharge,syscharge,up_volume,low_volume,atoms,probe,p30,p31,water = parse_space( pqrlines, potlines, box_array, center, nlines, E30_resid, K31_resid )

    max_volume=product(box_array)
    
    fieldfile = open( potentialfile.replace(".txt",".field"), "w")
    fieldfile.write("# %s\n" %potentialfile)
    fieldfile.write("#%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n"\
                    %("Linear fit","Poly fit", "Poly infl", "Sum Charge", "Sum +Charge", "Sum -Charge", "Sys Charge", "Max Volume", "Min Volume", "% Max Volume", "% Min Volume","Total Atoms", "% Atoms"))
    fieldfile.write(" %12.8f %12.8f %12.8f %12.3f %12.3f %12.3f %12.3f %12.4f %12.4f %12.4f %12.4f %12i %12.4f\n"\
                    %(linfield,polyfield,inflfield,charge,pcharge,ncharge,syscharge,up_volume,low_volume,up_volume/max_volume*100.,low_volume/max_volume*100.,atoms,atoms/(nlines-11.)*100.))
    
    write_array(probe,probe[0],fieldfile)
    if len(p30) != 0 :
        write_array(p30,p30[0],fieldfile)
    if len(p31) != 0 :
        write_array(p31,p31[0],fieldfile)
    if len(water) != 0 :
        write_array(water,water[0],fieldfile)
    fieldfile.close()

def parse_space( pqrlines, potlines, box, center, nlines, E30_resid, K31_resid ) :
    uppervolume=[]
    lowervolume=[]
    charge=[]
    pcharge=[]
    ncharge=[]
    syscharge=[]
    pprobe=[]
    p30=[]
    p31=[]
    water=[]
    mincoords,maxcoords=box_coords(center,box)
    
    for i in arange(nlines) :
        potential = float(potlines[i])
        line=pqrlines[i].split()
        partial=float(line[8])
        radius=float(line[9])
        syscharge.append(partial)
        if potential != 0 and "DUM" not in line :
            charge.append(partial)
            if partial >= 0 : pcharge.append(partial)
            elif partial < 0 : ncharge.append(partial)
            uppervolume.append(4./3.*pi*radius**2)
            lowervolume.append(atom_volume( radius, line, mincoords, maxcoords ))
            check_residue( line, "CNC", pprobe ) # nitrile probe resname, resid # changes
            check_residue( line, E30_resid, p30 ) 
            check_residue( line, K31_resid, p31 ) 
            check_residue( line, "SOL", water ) # the water inside the box

    return sum(charge), sum(pcharge), sum(ncharge), sum(syscharge), sum(uppervolume), sum(lowervolume), len(charge), pprobe, p30, p31, water

def atom_volume( radius, line, mincoords, maxcoords ) :
    coords=array( [float( line[i] ) for i in (5,6,7) ] )
    V = 4./3.*pi*radius**2
    Vt= V
    n=0
    for i in range(3) :
        hmin = abs( coords[i] - mincoords[i] )
        hmax = abs( coords[i] - maxcoords[i] )
        if hmin < radius :
            Vt -= (radius*hmin**2*pi - 1./3.*hmin**3*pi)
            n += 1
        elif hmax < radius :
            Vt -= (radius*hmax**2*pi - 1./3.*hmax**3*pi)
            n += 1
    if n == 0 or n == 1 : return Vt
    elif n == 2 and Vt < V/4. : return V/4.
    elif n == 3 and Vt < V/8. : return V/8.
    else : return Vt

    
def check_residue( pqrline, residue, an_array ) :
    try :
        ( residue ) 
        if "%s" %int(residue) == pqrline[4] :
            if len(an_array) == 0 : an_array.append(pqrline[3])
            an_array.append("%s.%s" %(pqrline[2],pqrline[4]))
    except :
        if residue == pqrline[3] :
            if len(an_array) == 0 : an_array.append(pqrline[3])
            an_array.append("%s.%s" %(pqrline[2],pqrline[4]))

def write_array( an_array, header, writefile ) :
    writefile.write("# %7s (%2i): "%(an_array[0], len(an_array)-1))
    for value in an_array[1:] :
        if value == an_array[-1] : writefile.write("%10s\n"%value)
        else : writefile.write("%10s "%value)

def coordComp(a,b):
	same=True
	for i in xrange(len(a)) :
		if a[i] != b[i] :
			return False
		else :
			same=True
	return True
def atomname2gmx( topfile ) :
    try :
        top=open(topfile)
    except :
        print "Coult not open topology file %s, exiting" %topfile
        sys.exit()
    name2gmx={}
    for line in top.readlines() :
        try :
            resname=line.split()[3]
            atom=line.split()[4]
            mass=float(line.split()[7])

            name2gmx[ resname, atom ] = mass
        except IndexError : pass # means this is not a line for an atom
        except ValueError : pass # means column 7 is not a mass
    top.close()
    return name2gmx

def entire_protein_box( pqrfile ):
    try :
        pqr=open(pqrfile)
    except :
        print "Could not open pqr file %s, exiting" %pqrfile
        sys.exit()
    pqrlines=pqr.readlines()
    pqr.close()

    xs=[]
    ys=[]
    zs=[]
    rs=[]
    for line in pqrlines :
        l=line.split()
        xs.append(float(l[5]))
        ys.append(float(l[6]))
        zs.append(float(l[7]))
        rs.append(float(l[9]))
    r=max(rs)
    minx=floor(min(xs)-r)
    maxx=ceil(max(xs)+r)
    miny=floor(min(ys)-r)
    maxy=ceil(max(ys)+r)
    minz=floor(min(zs)-r)
    maxz=ceil(max(zs)+r)

    xd=(maxx-minx)
    yd=(maxy-miny)
    zd=(maxz-minz)
    boxdime=array((xd,yd,zd))
    center=array(((maxx+minx)/2,(maxy+miny)/2,(maxz+minz)/2))

    return center, boxdime

def scn_box( pqrfile ) :
    try :
        pqr=open(pqrfile)
    except :
        print "Could not open pqr file %s, exiting" %pqrfile
        sys.exit()
    pqrlines=pqr.readlines()
    pqr.close()

    xs=[]
    ys=[]
    zs=[]
    rs=[]
    for line in pqrlines :
        if "CNC" in line :
            if "SG" in line or "CD" in line or "NE" in line :
                l=line.split()
                xs.append(float(l[5]))
                ys.append(float(l[6]))
                zs.append(float(l[7]))
                rs.append(float(l[9]))

    r=max(rs)
    minx=floor(min(xs)-r)
    maxx=ceil(max(xs)+r)
    miny=floor(min(ys)-r)
    maxy=ceil(max(ys)+r)
    minz=floor(min(zs)-r)
    maxz=ceil(max(zs)+r)

    xd=(maxx-minx)
    yd=(maxy-miny)
    zd=(maxz-minz)
    boxdime=array((xd,yd,zd))
    center=array(((maxx+minx)/2,(maxy+miny)/2,(maxz+minz)/2))

    return center, boxdime
