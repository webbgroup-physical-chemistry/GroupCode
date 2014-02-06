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
projects = (4, )
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
            #numeric_frame = int( framename.replace("frame","" ) )

            snapshot = 0
            for line in data :
                if line[0] != "@" and line[0] !="#" :
                    time, angle = line.split()
                    diff = abs( float(angle) - targetangle )        
                    if diff < smallest_diff and snapshot != 193 and snapshot != 178 :
                        smallest_diff = diff
                        closest_to_target = ( run, framename, snapshot, time, angle )
                    snapshot += 1

        print "\tfor", targetangle, "degrees, closest was run=%d framename=%s snapshot=%d time=%s angle=%s" % closest_to_target

        # int, string, int, string
        run, framename, snapshot, time = closest_to_target[:4]
        framenum = int( framename.replace( "frame", "" ) )

        # the snapshot is relative to the trajectory segment it's in, but we need the 
        # frame number in the whole trajectory
        # monomer:
        #   runs 0-8:
        #       frame0 = 201 snapshots (trivial case)
        #   runs 9-11: 
        #       frame1 = 201 snapshots
        #       frame2 = 201 snapshots
        #       frame3 = 201 snapshots
        # complex: 
        #   frame1 = 101 snapshots
        #   frame2 = 101 snapshots
        #   frame3 = 201 snapshots 
        if project < 11 : # monomer
            if run < 9 : # 201 snapshots in one frame, we're good to go
                pass
            else : # 201 frames/trajectory
                snapshot += (framenum-1)*201
        else : # complex
            snapshot += (framenum-1)*101

        # ../electrostatics-1b/project1/run0/frame165_C.fld.txt
        prdir = "../electrostatics-1b/project%d/run%d" % ( project, run )
        cfldfile = "%s/frame%d_C.fld.txt" % ( prdir, snapshot )
        gfldfile = "%s/frame%d_G.fld.txt" % ( prdir, snapshot )
        
        coutname = "fieldfiles/p%d_frame%d_C.fld.txt" % ( project, snapshot )
        goutname = "fieldfiles/p%d_frame%d_G.fld.txt" % ( project, snapshot )
        ccmd = "cp %s %s" % ( cfldfile, coutname )
        gcmd = "cp %s %s" % ( gfldfile, goutname )

        for cmd in ccmd, gcmd : 
            print cmd
            system(cmd)
