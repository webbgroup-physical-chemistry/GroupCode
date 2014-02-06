#!/usr/local/bin/python

"""
1. Read radii from AMBER.DAT.
    - all carbons are the same size
    - all nitrogens are the same size
    - etc.
    => store hydrogens by atom type
2. Read ffamber03.atp to translate GMX atom type "amber99_XX" to AMB atom type
3. Read topology file to (1) get charges and (2) translate GMX atom types to 
    residue+atom names
4. Read the pdb file and insert <charge> and <radius> columns
5. If use_solvent, for 'SOL' atoms, put tip3p charges and for sodium put +1


"""
use_solvent = True
tip3p_ocharge = -0.834
tip3p_hcharge = 0.417
sodium_charge= 1.000

import sys

try :
    dat, atp, top, pdb = sys.argv[1:]
except :
    print "usage: %s <dat> <atp> <top> <pdb>" % sys.argv[0] 
    sys.exit()

def readfile(filename):
    handle = open(filename)
    lines = handle.readlines()
    handle.close()
    return lines

datlines = readfile(dat)
atplines = readfile(atp)
toplines = readfile(top)
pdblines = readfile(pdb)

# 1. using the dat file, build a map from AMBER atom type to size
atomtype2size={}
for line in datlines :
    try :
        resname, atomname, charge, radius, atomtype = line.split()
        radius = float(radius)

        atomtype2size[ atomtype ] = radius # all atom type radii are the same
    except ValueError : pass # in case the line is not 5 things 

# 2. build a GMX atom type to AMB atom type dictionary
# this is very hackish
gmx2amb = {}
for line in atplines :
    parts = line.split()
    try :
        gmxname = parts[0]
        ambname = parts[3]
        if gmxname[-1] == "Y" : ambname = gmxname[-2:] # for CY and NY
        if ambname == "H0" : ambname = "H1" # ??

        gmx2amb[ gmxname ] = ambname
    except IndexError : pass

#for key in gmx2amb.keys() :
#    print key, gmx2amb[key]

# 3. read topology to get charges and build the atom name to GMX atom type map
atomname2gmx = {}
charges = {}
for line in toplines :
    try:
        parts=line.split()
        atomtype = parts[1]
        resname = parts[3]
        atomname = parts[4]
        charge = float( parts[6] )
    
        atomname2gmx[ resname, atomname ] = atomtype
        charges[ resname, atomname ] = charge
    except IndexError: pass # means this is not a line for an atom
    except ValueError: pass # means col 6 is not a charge

#for key in atomname2gmx.keys():
#    print key, atomname2gmx[key]
#for key in charges.keys():
#    print key, charges[key]

# 4. read the pdb and insert charges and radii
for line in pdblines :
    try :
        parts = line.split()
        atomname = parts[2]
        resname = parts[3]

        # crude, but the top has, eg, HD11 whereas the pdb has 1HD1
        if len(atomname)==4 :
            atomname = atomname[1:] + atomname[0]

        try:
            gmxname = atomname2gmx[ resname, atomname ]
            ambname = gmx2amb[ gmxname ]

            if ambname == "CY": ambname = "CT"  # this is a hack, but all C have the same size
            if ambname == "NY": ambname = "N" # another hack, same reason
            radius = "%8.3f" % atomtype2size[ ambname ] 

            charge = "% 8.4f" % charges[ resname, atomname ]
        
            print line[:54] + charge + radius
        except KeyError :
            radius = "%8.3f" % 0.0
            if resname == "Na+" and use_solvent :
                charge = "% 8.4f" % 1.0
                print line[:54] + charge + radius
            elif resname == "SOL" and use_solvent :
                if atomname == "OW" :
                    charge = "% 8.4f" % tip3p_ocharge
                    radius ="%8.3f" % 1.6612
                else : # hydrogen
                    charge = "% 8.4f" % tip3p_hcharge 
                print line[:54] + charge + radius
            else :
                pass
            
                #print line[:54], "could not translate this atom"
                #raise

    except IndexError : pass # not an atom line
