#!/usr/bin/python

# For each mutant (project) return the snapshot with CB-SG torsion closest 
# to the ideal alkane torsions of -180, -60, and 60 degrees.

from glob import glob
import sys
from os import system
from os.path import basename
from shutil import copy

targets = ( -180.0, -60.0, 60.0 )
projects = range(22)

for project in projects  :

    # get xvg files
    xvgfiles = glob( "project%d/run*/frame*xvg" % project ) 
    print "found %d xvg files for project %d" % ( len(xvgfiles), project)

    for targetangle in targets :
        smallest_diff = 100 # initialize
        for xvgfile in xvgfiles :
            file = open( xvgfile )
            data = file.readlines()
            file.close()

            ppart, rpart, fpart = xvgfile.split("/")
            run = int( rpart.replace("run", "" ) )
            framename = fpart.split(".")[0] 
            numeric_frame = int( framename.replace("frame","" ) )
    
            frame = 0
            for line in data :
                if line[0] != "@" and line[0] !="#" :
                    time, angle = line.split()
                    diff = abs( float(angle) - targetangle )        
                    if diff < smallest_diff :
                        smallest_diff = diff
                        closest_to_target = ( run, framename, frame, angle )
                    frame += 1

        print "\tfor", targetangle, "degrees, closest was run=%d framename=%s framenum=%d time=%s angle=%s" % closest_to_target
        
        # int, string, int, string
        run, framename, framenum, time = closest_to_target[:4]

        # ../electrostatics-1b/project1/run0/frame165_C.fld.txt
        prdir = "../project%d/run%d" % ( project, run )
        xtc = "%s/%s.xtc" % ( prdir, framename )
        tpr = "%s/frame0.tpr" % prdir 

        outname = "grofiles/p%d_%d.gro" % ( project, targetangle )
        logname = "log_p%d_%d.txt" % ( project, targetangle )
        trjconvcmd = "echo 0 | trjconv -f %s -s %s -dump %s -o %s >& %s" % ( xtc, tpr, time, outname, logname )
        print trjconvcmd
        system(trjconvcmd)
