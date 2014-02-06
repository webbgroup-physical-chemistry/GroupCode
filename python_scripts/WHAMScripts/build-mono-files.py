#!/usr/bin/python
#../simready-analysis/G28C-0.xvg.txt    0   45  1000 0.400907
# runs 0-8   frame0.angle.xvg
# runs 9-11  frame1.angle.xvg
#            frame2.angle.xvg
#            frame3.angle.xvg

mutantlist=open("../mutantlist.txt").readlines()

dphi = 45
kfac = 1000
beta = .400907

# make a map for runs 9-11
file = open( "mutant_anglelist" )
pr2angle = {}
for line in file.readlines() :
    p, r, angle = line.split()
    p = int(p)
    r = int(r)
    angle = float(angle)
    pr2angle[ p, r ] = angle

for project in range(11) :
    mutant = mutantlist[project].strip()
    filename = "mono.%s.files" % mutant
    file = open( filename, "w" )

    for run in range(8) :
        xvg="project%d/run%d/frame0.angle.xvg" % ( project, run )
        angle = 45 * run
        file.write( "%25s\t%d\t%d\t%4.3f\t%6.6f\n" % ( xvg, angle, dphi, kfac, beta ) )

    xvg = "project%d/run8/frame0.angle.xvg" % project
    angle = 0 
    file.write( "%25s\t%d\t%d\t%4.3f\t%6.6f\n" % ( xvg, angle, dphi, 0, beta ) )

    for run in range(9,12) :
        angle = pr2angle[ project, run ]
        for xvgnumber in range(1,4) :
            xvg="project%d/run%d/frame%d.angle.xvg" % ( project, run, xvgnumber )
            file.write( "%25s\t% 5d\t%d\t%4.3f\t%6.6f\n" % ( xvg, angle, dphi, kfac, beta ) )
            
        
