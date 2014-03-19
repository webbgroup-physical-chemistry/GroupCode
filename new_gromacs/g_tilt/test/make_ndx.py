#! /anaconda/bin/python

import sys
sys.path.append("/Users/ritchie/Utilities/python_scripts")
from structure_manip import *

"""
LFD=structure("1LFD.gro")
#gro=structure("Rap_E30_K31+N27C.60-300.18.gro")
gro=structure("nosol.gro")

s1,s2,n1,n2=align_sequence(LFD,gro)
print n1

LFD_Ral = []
LFD_GTPase = []
gro_Ral = []
gro_GTPase = []
for i in range(len(LFD.index)) : 
    if LFD.resname[i] != "HOH" :
        if LFD.atom[i] == "CA" : 
            if LFD.resid[i] <= 100 : 
                LFD_Ral.append(LFD.index[i])
            else :
                LFD_GTPase.append(LFD.index[i])

for i in range(len(gro.index)) : 
    if gro.resname[i] != "SOL" and gro.resname[i] != "Na+" : 
        if gro.atom[i] == "CA" :
            if gro.resid[i] <= 100 : 
                gro_Ral.append(gro.index[i])
            else :
                gro_GTPase.append(gro.index[i])
"""
new=structure("Rap_alignd_to_1lfd.pdb")
for i in range(len(new.line)) :
    print new.line[i].strip()
    print new.resid[i]
