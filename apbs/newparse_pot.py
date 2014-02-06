#! /usr/bin/python

import sys
from numpy import *
from numpy.linalg import norm
from optparse import OptionParser
from pylab import *

parser=OptionParser()
parser.add_option("-o", "--output", dest="out", help="Output file name", default="out.pot")
parser.add_option("-f", "--focus", dest="strategy", help="Focusing strategy", default="G")
parser.add_option("-d", "--dimension", dest="dime", help="Dimension", default="193")
parser.add_option("-n", "--nframes", dest="nframes", help="Number of frames", default="600")
parser.add_option("-r", "--order", dest="order", help="Order of the polynomial fit", default="10")
parser.add_option("-p", "--position", dest="position", help="Fraction of bond length to compute the field at, default: .5", default=".5")
parser.add_option("-a", "--ndummy",dest="ndummy", help="Number of dummy atoms",default="11")
(option,args) = parser.parse_args()
                
out=option.out
strategy=option.strategy
dime=int(option.dime)
nframes=int(option.nframes)+1
order=int(option.order)
position=float(option.position)
ndummy=int(option.ndummy)

def getBondLength( structurefile ) :
    for line in structurefile.readlines() :
        if "CNC" in line and "NE" in line :
            nx = float(line[32:38])
            ny = float(line[40:46])
            nz = float(line[48:54])
        if "CNC" in line and "CD" in line :
            cx = float(line[32:38])
            cy = float(line[40:46])
            cz = float(line[48:54])
    BondVector = ((nx-cx, ny-cy, nz-cz))
    return norm(BondVector)

def getPot( structurefile, potentialfile, order, position) :
    bl = getBondLength( structurefile )
    iteration = 0
    pot = []
    pos = []
    start=(ndummy-1)/2-2
    end=(ndummy-1)/2+3
    print "Looking at potentials from atom %i to atom %i." %(start,end)
    for potvalue in potentialfile.readlines()[4:4+ndummy] :
        dist = bl*iteration/(ndummy-1)
        potential = float(potvalue.split()[0])
        if potential == 0 :
            error = True
        pot.append(potential)
        pos.append(dist)
        iteration +=1
    if order == 1 :
        polynomial=poly1d(polyfit(pos[start:end],pot[start:end],order))
    else :
        polynomial = poly1d(polyfit(pos,pot,order))      
    deriv = polyder(polynomial)
    final = -deriv(position*bl)
    seq=arange(.4*bl,.65*bl,.001)
    inflection_point =  min(seq, key=lambda x: -deriv(x))
    r_distance = inflection_point/bl
    inflection_field = min(map(lambda x: -deriv(x),seq))
    #if order == 1 :
    #    print -deriv(pos[5])
    #else :
    #    for xs in pos[3:8] :
    #        print -deriv(xs)
    #print "\n"
    #xp=linspace(0,pos[10],100)
    #plot(pos, pot, "xr", xp, polynomial(xp), "-b")
    return final, inflection_field, inflection_point, r_distance
    
def parseFiles( n, focus, dime, order, position, out ) :
    potential = []
    inflection_potential = []
    for i in range(n) :
        pqrname = "frame%d.pqr" %i
        apbsname = "frame%d_%s_%d.txt" %( i, focus, dime )
        try :
            pqrin = open( pqrname )
            try :
                apbsin = open( apbsname )
            except :
                print "*****Missing %s, skipping...*****" %apbsname
                pass
            frame_potential, inflection_field, inflection_point, r_dist = \
                             getPot(pqrin, apbsin, order, position)
            if frame_potential == 0 :
                frame_potential = "Error"
            potential.append(frame_potential)
            inflection_potential.append(inflection_field)
            pqrin.close()
            apbsin.close()
        except :
            print "*****Missing %s, skipping...*****" %pqrname
            potential.append("Error %s" %i )
            inflection_potential.append("Error %s" %i )
            pass
    saveOutput(potential, out)
    saveOutput(inflection_potential, "inflection_%s" % out)
    return potential, inflection_potential

def saveOutput( AnArray, out ) :
    outfile = open( out, "w" )
    for element in AnArray :
        outfile.write( "%s\n" %element )
    outfile.close()
        

potential, inflection = parseFiles( nframes, strategy, dime, order, position, out)



