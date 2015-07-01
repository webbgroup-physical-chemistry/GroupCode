#! /usr/bin/env python

import numpy
import os
from optparse import OptionParser
import textwrap
import sys

""" 
template.in 
"""
template = "\
\n\
read \n\
    mol pqr PQRNAME \n\
end \n\
\n\
# Solvation Energy - SOLVATED STATE \n\
elec name solv \n\
    mg-manual \n\
    dime DIME \n\
    glen GLEN \n\
    gcent mol 1 \n\
    mol 1 \n\
    lpbe \n\
    bcfl mdh \n\
    ion charge 1 conc CONC radius 2.0 \n\
    ion charge -1 conc CONC radius 2.0 \n\
    pdie PDIE \n\
    sdie SDIE \n\
    chgm spl2 \n\
    srfm mol \n\
    srad 1.4 \n\
    sdens 10.0 \n\
    swin 0.3 \n\
    temp TEMP \n\
    calcenergy total \n\
    calcforce no \n\
    write atompot flat out_solvated \n\
end \n\
\n\
# Solvation Energy - REFERENCE STATE \n\
elec name ref  \n\
    mg-manual \n\
    dime DIME \n\
    glen GLEN \n\
    gcent mol 1 \n\
    mol 1 \n\
    lpbe \n\
    bcfl mdh \n\
    ion charge 1 conc CONC radius 2.0 \n\
    ion charge -1 conc CONC radius 2.0 \n\
    pdie PDIE \n\
    sdie PDIE \n\
    chgm spl2 \n\
    srfm mol \n\
    srad 1.4 \n\
    sdens 10.0 \n\
    swin 0.3 \n\
    temp TEMP \n\
    calcenergy total \n\
    calcforce no \n\
    write atompot flat out_vacuum \n\
end \n\
\n\
quit\n\
"

"""
An option parser
"""
parser = OptionParser()
parser.add_option( "-f"\
                  ,"--pqrFile"\
                  ,dest = "pqrfile"\
                  ,help = "<.pqr> File name.  Do NOT use a .pdb!  This will search for space-delimited values, which may be inaccurate for a .pdb file"\
                  )
parser.add_option( "-i"\
                  ,"--inputFile"\
                  ,dest = "inputfile"\
                  ,help = "Name of the APBS input .in file"\
                  ,default = "input.in"\
                  )
parser.add_option( "-g"\
                  ,"--gridSpacing"\
                  ,dest = "gridTarget"\
                  ,help = "Space (Angstroms) between grid points in each dimension.  DEFAULT = 0.25"\
                  ,default = 0.25\
                  ,type="float"\
                  )
parser.add_option( "-c"\
                  ,"--ionConcentration"\
                  ,dest = "ionConc"\
                  ,help = "Ion concentration.  DEFAULT = 0.0"\
                  ,default = 0.0\
                  ,type="float"\
                  )
parser.add_option( "-p"\
                  ,"--soluteDielectric"\
                  ,dest = "pdie"\
                  ,help = "Solute (protein) dielectric.  DEFAULT = 1.0"\
                  ,default = 1.0\
                  ,type="float"\
                  )
parser.add_option( "-s"\
                  ,"--solventDielectric"\
                  ,dest = "sdie"\
                  ,help = "Solvet dielectric.  DEFAULT = 78.0"\
                  ,default = 78.0\
                  ,type="float"\
                  )   
parser.add_option( "-t"\
                  ,"--temperature"\
                  ,dest = "temp"\
                  ,help = "Temperature, in Kelvin.  DEFAULT = 300.0"\
                  ,default = 300.0\
                  ,type="float"\
                  )       
parser.add_option( "-S"\
                  ,"--StripPath"\
                  ,dest = "StripPath"\
                  ,help = "Strip the pathname from the pqr file name in input.in. DEFAULT = True"\
                  ,action = "store_false"\
                  ,default = True\
                  )
(options,args) = parser.parse_args()
pqrfile = options.pqrfile
infile = options.inputfile
gridTarget = options.gridTarget
ionConc = options.ionConc
pdie = options.pdie
sdie = options.sdie
temp = options.temp
StripPath = options.StripPath


"""
Adding some general file utilities I like
"""
maxwidth = 65

def readFile( filename, kill=True ) :
    try :
        File = open(filename)
        FileLines = File.readlines()
        File.close()
    except :
        print "Cannot open %s"%filename
        if kill :
            print "EXITING..."
            sys.exit()
        return 0
    return FileLines

def printbox( string ) :
    printw( "\n"+"-"*maxwidth )
    printw("|"+" "*(maxwidth-2)+"|")
    for line in textwrap.wrap(string,width=maxwidth-2) :
        printw("|"+line.center(maxwidth-2)+"|")
    printw("|"+" "*(maxwidth-2)+"|")
    printw( "-"*maxwidth +"\n")
    return True

def printcenter( string ) :
    for line in textwrap.wrap(string,width=maxwidth) :
        printw(line.center(maxwidth))
    return True

def printw( string ) :
    lines = string.split('\n')
    for line in lines :
        print textwrap.fill(line,width=maxwidth)
    return True
    
def backup_outname( filename ) :
    filename = os.path.abspath(filename)
    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    copyname = filename
    n = 1
    if os.path.isfile(filename) :
        while os.path.isfile(copyname) :
            if n == 100 :
                printw( "Will no make 100 copies of %s, exiting..."%filename )
                return
            copyname = "%s/#%s.%i#"%(dirname,basename,n)
            n += 1
        os.rename(filename, copyname)
        printbox( "Backing up %s to %s"%(filename,copyname))
    return filename

def writeLines( name, newlines ) :
    if len(name) > 1054 :
        printw("The name is REALLY long!  You probably put the line string there by mistake.  Edit the code if you want this to work...")
        return 0
    backup_outname( name )
    File = open(name,"w")
    File.write(newlines)
    File.close()
    printw("\n>> Done Writing %s <<\n"%name)
    return name

def GetDime( glen ) :
    dime = []
    for j in range(3) :
        each = glen[j]
        d = each/gridTarget + 1
        if d < 33. : 
            printbox("There are fewer than 33 grid points in dimension %i; setting to 33.  You may consider decreasing the grid point spacing"%(j+1))
            d = 33.
        else : 
            for i in range(2,200) :
                if d < i*2**5+1 : 
                    d = int((i-1)*2**5+1)
                    break
        dime.append(d)
    return dime

def GenerateDimensions():
    xs = []
    ys = []
    zs = []
    radii = []
    pqrLines = readFile(pqrfile)
    for line in pqrLines : 
        if line.startswith("ATOM") or line.startswith("HETATM") : 
            coords = line.split()[5:10]
            xs.append(float(coords[0]))
            ys.append(float(coords[1]))
            zs.append(float(coords[2]))
            radii.append(float(coords[4]))
    dimX = max(xs) - min(xs)
    dimY = max(ys) - min(ys)
    dimZ = max(zs) - min(zs)
    maxRadii = max(radii)
    glen = numpy.ceil(max([ dimX,dimY,dimZ ]) + 2*maxRadii)
    dime = GetDime([glen,glen,glen])    # Making it a cubic box.
    return ([glen,glen,glen],dime)

def writeInput(glen,dime) : 
    newInput = template.replace("CONC","%.4f"%ionConc)
    newInput = newInput.replace("DIME","%i %i %i"%(dime[0],dime[1],dime[2]))
    newInput = newInput.replace("GLEN","%.3f %.3f %.3f"%(glen[0],glen[1],glen[2]))
    newInput = newInput.replace("SDIE","%.3f"%(sdie))
    newInput = newInput.replace("PDIE","%.3f"%(pdie))
    newInput = newInput.replace("TEMP","%.3f"%(temp))
    pqrname = pqrfile
    if StripPath : 
        pqrname = os.path.basename(pqrfile)
    newInput = newInput.replace("PQRNAME","%s"%pqrname)
    writeLines(infile,newInput)            
        
if __name__ == "__main__":
    (glen,dime) = GenerateDimensions()
    writeInput(glen,dime)
    #writeLines("/Users/ritchie/Utilities/apbs/solv.template.in",template)
    
    
