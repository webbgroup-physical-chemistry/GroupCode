#! /usr/local/bin/python

atp="/usr/local/gromacs/share/gromacs/top/ffamber03a.atp"
f=open(atp)
fl=f.readlines()
f.close()

mapping={}
for line in fl :
    if not line.startswith(";") :
        l=line.split()
        if "amber99_" in l[0] :
            mapping[l[0]]=l[3]

fff="/Users/ritchie/software/gromacs-4.5.5/share/top/amber03a.ff/aminoacids.rtp"
ff=open(fff)
fl=ff.readlines()
ff.close()


new=open(fff.replace(".rtp","_a.rtp"),"w")

for line in fl :
    for m in mapping :
        line=line.replace(" %s "%m," %s "%mapping[m])
    new.write(line)
new.close()

