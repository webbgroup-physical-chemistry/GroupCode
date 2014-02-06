import sys

try :
    filename=sys.argv[1]
    file=open(filename)
    filelines=file.readlines()
    file.close
except :
    print "usage: %s pqr" %sys.argv[0]
    sys.exit()

water=""
protein=""
for line in filelines :
    if "CNC" in line :
        if "CD" in line or "NE" in line or "SG" in line :
            line=line.replace(line[55:62],"%7.4f" %0)
    water+=line
    if "SOL" not in line :
        protein+=line

watername=filename.replace(".pqr","_nc_w.pqr")
proteinname=filename.replace(".pqr","_nc.pqr")

f1=open(watername,"w")
f1.write(water)
f1.close()

f2=open(proteinname,"w")
f2.write(protein)
f2.close()
