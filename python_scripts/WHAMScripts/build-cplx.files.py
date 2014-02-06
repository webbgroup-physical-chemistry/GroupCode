#!/usr/bin/python
# $ ls ../project11/run0/
# frame1.angle.xvg frame2.angle.xvg frame3.angle.xvg

mutantlist=open("../mutantlist.txt").readlines()

dphi = 45
kfac = 1000
beta = .400907

for project in range(11,22) :
    mutant = mutantlist[project-11].strip()
    filename = "cplx.%s.files" % mutant
    file = open( filename, "w" )

    for run in range(6) :
        for xvgnumber in (1,2,3):
            xvg="project%d/run%d/frame%d.angle.xvg" % ( project, run, xvgnumber )
            angle = 60 * run
            file.write( "%25s\t%d\t%d\t%4.3f\t%6.6f\n" % ( xvg, angle, dphi, kfac, beta ) )
