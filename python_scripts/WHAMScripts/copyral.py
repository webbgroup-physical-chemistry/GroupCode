#!/usr/bin/python
# this copies files for projects 0-10 , runs 9-11 to the appropriate directory and
# renames

from shutil import move
file=open( "mutant_anglelist" )
filedata=file.readlines()
file.close()

for line in filedata :
    dir, ma = line.split()

    for i in (1,2,3) :
        move( "%s.%d.xvg" % ( ma, i ), "%s/frame%d.angle.xvg" % ( dir, i ) )             

