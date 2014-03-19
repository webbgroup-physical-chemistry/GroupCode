#! /anaconda/bin/python


gro = "Rap_alignd_to_1lfd.pdb"
readgro = open(gro)
grolines = readgro.readlines()
readgro.close()

DangerZone = False

for line in grolines :
    if "ATOM" in line :
        resid = line[22:27].split()[0]
        if DangerZone and "GNP" not in line : 
            print line.replace(line[22:27],"%5i"%newresid).strip()
            if "Na+" in line or "HW2" in line : 
                newresid += 1
        if not DangerZone or (DangerZone and "GNP" in line) : 
            print line.strip()
        if not DangerZone and int(resid) == 272 :
            DangerZone = True
            newresid = 273
print "0.00 0.00 0.00"
