#! /usr/bin/python

from numpy import array, pi, cos, sin
from optparse import OptionParser
from math import atan2
import os
import sys

parser = OptionParser()
parser.add_option("-A", "--azimxvg", dest="azimxvg", help="<azimuthal>.xvg \
file containing data per frame")
parser.add_option("-P", "--polarxvg", dest="polarxvg", help="<polar>.xvg \
file containing data per frame")
parser.add_option("-a", "--azimuthal", dest="azim", help="boltzmann weighted \
azimuthal angle output file")
parser.add_option("-p", "--polar", dest="polar", help="boltzmann weighted \
polar angle output file")
parser.add_option("-Z", "--scazimxvg", dest="scazimxvg", help="<side chain \
azimuthal>.xvg file containing data per frame", default="")
parser.add_option("-z", "--scazimuthal", dest="scazim", help="boltzmann weighted \
side chain azimuthal angle output file", default="")
parser.add_option("-L", "--scpolarxvg", dest="scpolarxvg", help="<side chain \
polar>.xvg file containing data per frame", default="")
parser.add_option("-l", "--scpolar", dest="scpolar", help="boltzmann weighted \
side chain polar angle output file", default="")
parser.add_option("-f", "--xtc", dest="xtc", help="Trajectory: xtc trr",\
                  default="")
parser.add_option("-s", "--gro", dest="gro", help="Structure+mass: gro pdb",\
                  default="")
parser.add_option("-o", "--out", dest="outname", help="<out>.pdb", \
                  default="out")
parser.add_option("-e", "--verror", dest="error", help="Factor to increase \
allowed side chain vector variance by.  Default = 2", default=2)
parser.add_option("-v", "--perror", dest="perror", help="Factor to increase \
allowed primary vector variance by.  Default = 1", default=1)
options, args = parser.parse_args()

azimxvg = options.azimxvg
polarxvg = options.polarxvg
azim = options.azim
polar = options.polar
outname = options.outname
scazimxvg = options.scazimxvg
scazim = options.scazim
scpolarxvg = options.scpolarxvg
scpolar = options.scpolar
xtc = options.xtc
gro = options.gro
error_range = float(options.error)
n_error_range = float(options.perror)



if len(xtc) > 0 and len(gro) > 0 :
    make_pdbs = True
    print "If a representative structure is found, a pdb file will be \
generated\n"
else :
    make_pdbs = False
    print "If a representative structure is found, a pdb file will NOT \
be generated\n"

if len(scpolarxvg)*len(scpolar)*len(scazim)*len(scazimxvg) > 0 :
    use_side_chain = True
    print "Looking for representative side chain structure as well...\n"
else :
    use_side_chain = False


def handle_angles( theta ) :
    x, y = cos(theta*pi/180), sin(theta*pi/180)
    return int( atan2(y,x)*180/pi )

mean = {}
std = {}

positions = []

for parameter in array((azim, polar, scazim, scpolar)) :
    if len(parameter) > 0 :
        positions.append(parameter)

for angle in array((positions)) : 
    try :
        line = open(angle).readlines()[-1]
        mean[angle] = handle_angles( float(line.split()[2]) )
        if angle == azim or angle == polar :
            std[angle] = handle_angles( float(line.split()[3]) )*n_error_range
        elif angle == scazim or angle == scpolar :
            std[angle] = handle_angles( float(line.split()[3]) )*error_range

    except :
        print "\n\n*****Error trying to open Boltzmann weighted output file, %s*****\n\n" %angle
        sys.exit()
if use_side_chain :
    print "Looking for an azimuthal angle of %i +/- %i and a polar angle of %i \
+/- %i.  The side chain azimuthal angle is %i +/- %i and the polar angle is %i \
+/- %i.\n" \
%(mean[azim], std[azim], \
  mean[polar], std[polar], \
  mean[scazim], std[scazim], \
  mean[scpolar], std[scpolar])
else :    
    print "Looking for an azimuthal angle of %i +/- %i and a polar angle of %i \
+/- %i.\n" %(mean[azim], std[azim], \
             mean[polar], std[polar])


frames = {}
for angle in array((positions)) :
    matches = 0

    if angle == azim :
        looking_at = azimxvg
    elif angle == polar :
        looking_at = polarxvg
    elif angle == scazim :
        looking_at = scazimxvg
    elif angle == scpolar :
        looking_at = scpolarxvg

    error = std[angle]           
    print "Parsing through %s" %looking_at
        
    try :
        file = open(looking_at)
        lines = file.readlines()
    
        for line in range( len(lines) ) :
            frame = handle_angles( float(( lines[line].strip() )) ) 
            if mean[angle] - error <= frame \
               <= mean[angle] + error :
                matches += 1
                key = "%i.%s" % (line, angle)
                frames[key] = ((frame, abs(frame-mean[angle])))

        else :
            print "Found %i matching angles..." % matches
    except :
        print "\n\n*****Error trying to open %s...*****\n\n" %looking_at
        sys.exit()


try :
    variances = {}
    var_sum = []

    for i in range(len(lines)) :
        try :
            akey = "%i.%s" %(i, azim)
            pkey = "%i.%s" %(i, polar)
            if use_side_chain :
                sakey = "%i.%s" %(i, scazim)
                spkey = "%i.%s" %(i, scpolar)
                aframe, pframe, saframe, spframe =\
                        frames[akey], frames[pkey], frames[sakey], frames[spkey]
            else :
                aframe, pframe = frames[akey], frames[pkey]
                saframe = spframe = ((0,0))
                
            variances[i] = aframe[0], pframe[0], saframe[0], spframe[0], \
                        (aframe[1]**2+pframe[1]**2\
                         +saframe[1]**2+spframe[1]**2)**.5
            var_sum.append((aframe[1]**2+pframe[1]**2+\
                            saframe[1]**2+spframe[1]**2)**.5)

        except :
            pass

    if len(var_sum) == 0 :
        print "\nThere are no representative frames in %s, %s.\n\n" %( azimxvg, \
                                                                  polarxvg )
        sys.exit()

    representative_frames = []
    for i in variances :
        if variances[i][4] == min(var_sum) :

            if use_side_chain :
                print "\nFrame %i is representative with an azimuthal angle" %i,
                print "of %i and a polar angle of %s." %(variances[i][0],\
                                                     variances[i][1]),
                print "The side chain azimuthal angle is %s and the polar \
is %s.  The root distance is %.2f. \n\n" \
%( variances[i][2], variances[i][3], variances[i][4] )
            else :
                print "\nFrame %i is representative with an azimuthal angle" %i+1,
                print "of %i and a polar angle of %s.The root distance is \
%.2f.\n\n" %(variances[i][0], variances[i][1], variances[i][4])
            representative_frames.append(i+1)

except :
    print "There are no representative frames in this trajectory"
    make_pdbs = False
    
if make_pdbs :
    make_ndx = "/home/ritchie/Utilities/NoSolIndex.py"
    framenumber = 0
    for frame in representative_frames :
        framenumber += 1
        time_1 = frame*5-1
        time_2 = frame*5+4
        cmd1 = "%s %s > NoSolvent.ndx" %( make_ndx, gro )
        cmd2 = "trjconv -f %s -s %s -o %s%s..pdb -n NoSolvent.ndx -sep -b %s -e %s -pbc nojump; rm NoSolvent.ndx"\
        %( xtc, gro, outname, framenumber, time_1, time_2)

        #print cmd1
        os.system(cmd1)
        #print cmd2
        os.system(cmd2)

#        cmd3 = "/home/ritchie/Utilities/ffamber_names.sh %s%s.0.pdb %s%s.pdb ;\
#rm %s%s.0.pdb "
#        os.system( cmd3 )
