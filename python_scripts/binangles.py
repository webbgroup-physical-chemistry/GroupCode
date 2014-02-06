#!/usr/bin/python

from numpy import arange, histogram
import sys

binwidth = 5 # degrees
binmin = -180
binmax = 180
bins = arange( binmin, binmax+binwidth, binwidth )

try:
    xvgfile = sys.argv[1]
    outname = sys.argv[2]
except:
    print "usage ./script <xvg> <outname>"
    sys.exit()

file = open(xvgfile)
xvglines=file.readlines()
file.close()

data = []
for line in xvglines :
    if line[0] != "#" and line[0] != "@" :
        time, angle = line.split()
        angle = float(angle)
        bin = int(( angle + 180 )/5 )
        if bin ==72: bin = 0
        data.append( "%d\n" % bin )
        #print time, angle, bin

file = open( outname, "w" )
file.writelines(data)
file.close() 
