#!/usr/bin/python
# this copies files for projects 11-21

from shutil import move
mutantlist = open("mutantlist.txt").readlines()

for project in range(11,22):
    mutant = mutantlist[project-11].strip()

    for run in range(6) :
        angle = run*60
        
        dir = "project%d/run%d" % ( project, run )
        xvg1 = "tmp/%s_%d.xvg" % ( mutant, angle )
        xvg2 = "tmp/%s_%d_2.xvg" % ( mutant, angle )
        xvg3 = "tmp/%s_%d_3.xvg" % ( mutant, angle )
   
        out1 = "%s/frame1.angle.xvg" % dir
        out2 = "%s/frame2.angle.xvg" % dir
        out3 = "%s/frame3.angle.xvg" % dir

        move( xvg1, out1 )
        move( xvg2, out2 )
        move( xvg3, out3 )
 
