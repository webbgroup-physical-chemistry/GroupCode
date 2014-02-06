#! /usr/bin/python
import sys
try: 
    file=open(sys.argv[1])
except:
    print "please provide a .xvg file from running g_angle -ov -type dihedral"

mybins=file.readlines()
for i in range(0,72):
    count=0
    frames=0
    for line in mybins:
        frames+=1.
        if int(line)==i:
            count += 1.
        if i==0:
            endcount=float(count/frames)        
    print i*5, float(count/frames) 
print "360",endcount
