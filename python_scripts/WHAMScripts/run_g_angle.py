#!/usr/bin/python
# run g_angle for projects 0-10, runs 0-8

#xtcname=../G28C-0_out.xtc
#ndxname=$mutant.ndx # need to make these from indexfiles.txt
#outputname=$mutant-$angle-avg.xvg

from os import system

mutantlist = open( "mutantlist.txt" ).readlines()

for project in range(11) :
    mutant = mutantlist[project].strip()
    ndxfile = "%s.ndx" % mutant
 
    for run in range(9) :
        dir = "project%d/run%d" % ( project, run )
        outfile = "%s/frame0.angle.xvg" % dir
        xtcfile = "../%s/frame0.xtc" % dir

        cmd = "g_angle -f %s -n %s -ov %s -noxvgr -type dihedral" % ( xtcfile, ndxfile, outfile )
        print cmd
        #system(cmd)
