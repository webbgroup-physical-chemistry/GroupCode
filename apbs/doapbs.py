#! /usr/local/bin/python

from apbs_setup import *
import os


from optparse import OptionParser
parser=OptionParser()
parser.add_option("-f",\
                  "--pqrfile",\
                  dest="pqrfile",\
                  help=".pqr file"\
                  )
parser.add_option("-t",\
                  "--top",\
                  dest="topfile",\
                  help=".top or .itp file (for mass information)"\
                  )
parser.add_option("-P",\
                  "--PATH",\
                  dest="PATH",\
                  help="Path to /bin/apbs"\
                  ,default="/usr/local/apbs/bin/apbs"\
                  )
parser.add_option("-n",\
                  "--norun",\
                  dest="runapbs",\
                  default=True,\
                  action="store_false",\
                  help="Do not run apbs, just analyze output"\
                  )
parser.add_option("-e",\
                  "--E30_resid",\
                  dest="E30",\
                  help="Residue numbers of position 30. Default=133 (Rap1a/Ral)",\
                  default="133"\
                  )
parser.add_option("-k",\
                  "--K31_resid",\
                  dest="K31",\
                  help="Residue numbers of position 31. Default=134 (Rap1a/Ral)",\
                  default="134"\
                  )
option,args=parser.parse_args()

pqrfile=option.pqrfile
topfile=option.topfile
PATH=option.PATH
runapbs=option.runapbs
E30_resid=float(option.E30)
K31_resid=float(option.K31)

name2gmx=atomname2gmx(topfile)

potfiles=[]
boxes=[]
centers=[]

pqr=open(pqrfile)
pqrlines=pqr.readlines()
nlines=len(pqrlines)
pqr.close()
SySCoM,CD,NE,BV=CoM(pqrlines, name2gmx)
BondLength=linalg.norm(BV)

protein_center, protein_box = entire_protein_box(pqrfile)
scn_center, scn_box = scn_box(pqrfile)
b19=array([19.,19.,19.])
b10=array([10.,10.,10.])
d161=array([161.,161.,161.])
d193=array([193.,193.,193.])

protein_dime = getdime(protein_box,19/96.)
scn_dime = getdime(scn_box,10/256.)

bcout,dxname=make_bc_input(pqrfile)

potfiles.append(make_apbs_input(pqrfile,dxname,CD,b19,d161,nlines,first=True,dowrite=runapbs)[0])
boxes.append(b19)
centers.append(CD)

potfiles.append(make_apbs_input(pqrfile,dxname,CD,b19,d193,nlines,dowrite=runapbs)[0])
boxes.append(b19)
centers.append(CD)

potfiles.append(make_apbs_input(pqrfile,dxname,CD,b10,d161,nlines,dowrite=runapbs)[0])
boxes.append(b10)
centers.append(CD)

if '_scn' in pqrfile :
    potfiles.append(make_apbs_input(pqrfile,dxname,scn_center,scn_box,scn_dime,nlines,dowrite=runapbs)[0])
    boxes.append(scn_box)
    centers.append(scn_center)
else :
    potfiles.append(make_apbs_input(pqrfile,dxname,protein_center,protein_box,protein_dime,nlines,dowrite=runapbs)[0])
    boxes.append(protein_box)
    centers.append(protein_center)

finalpot,final_inputfile=make_apbs_input(pqrfile,dxname,CD,b10,d193,nlines,dowrite=runapbs)
potfiles.append(finalpot)
boxes.append(b10)
centers.append(CD)

if runapbs :
    final_inputfile.write("quit\n")
    final_inputfile.close()
    try :
        bcoutdx=open(dxname)
        bcoutdx.close()
        print "%s already exists, skipping to second-stage calculation" %dxname
        cmd="%s %s" %(PATH,pqrfile.replace(".pqr",".in"))
    except :
        cmd="%s %s;%s %s" %(PATH,bcout,PATH,pqrfile.replace(".pqr",".in"))
    os.system(cmd)
else : print "Analyzing %s outputs." %pqrfile


for i in arange(len(potfiles)) :
    parse_pot(pqrlines, potfiles[i],BondLength,boxes[i],centers[i],nlines, E30_resid, K31_resid)    

print potfiles
print boxes
print centers
