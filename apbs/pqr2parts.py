#! /usr/local/bin/python

import sys
from numpy import array
from optparse import OptionParser
parser=OptionParser()
parser.add_option("-f",\
                  "--pqrfile",\
                  dest="pqrfile",\
                  default="12-E30D+N27C_180-180-frame25.pqr",\
                  help=".pqr file"\
                  )
options,args=parser.parse_args()
pqrfile=options.pqrfile

pqr=open(pqrfile)
pqrlines=pqr.readlines()
pqr.close()

max_ral=100
max_rap=272

def resid( line ) :
    return float(line[23:26].split()[0])
def zero_charge( line ) :
    return line.replace(line[55:62]," 0.0000")
def in_scn( line ) :
    if "CNC" in line :
        if "CD" in line or "NE" in line or "SG" in line :
            return True
    return False
def in_sol( line ) :
    if "SOL" in line :
        return True
    return False
def in_rap( line ) :
    if resid(line) > max_ral and not in_sol(line) :
        return True
    return False
def in_ral( line ) :
    if resid(line) <= max_ral and not in_sol(line) and not in_scn(line):
        return True
    return False

rap=open(pqrfile.replace(".pqr","_rap.pqr"),"w")
ral=open(pqrfile.replace(".pqr","_ral.pqr"),"w")
scn=open(pqrfile.replace(".pqr","_scn.pqr"),"w")
sol=open(pqrfile.replace(".pqr","_sol.pqr"),"w")

for line in pqrlines :
    if in_scn(line) : scn.write(line)
    else : scn.write(zero_charge(line))
    if in_sol(line) : sol.write(line)
    else: sol.write(zero_charge(line))
    if in_ral(line) : ral.write(line)
    else : ral.write(zero_charge(line))
    if in_rap(line) : rap.write(line)
    else : rap.write(zero_charge(line))


rap.close()
ral.close()
scn.close()
sol.close()
