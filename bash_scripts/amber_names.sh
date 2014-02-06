#! /bin/bash

pdb=$1

if [ -z $pdb ] ; then 
	echo "usage: $0 <pdb>"
	exit
else 
	echo Looking at $pdb
fi

sed -i .bu -e 's/NGL /NGLY/' $pdb
sed -i .bu -e 's/CTH /CTHR/' $pdb
sed -i .bu -e 's/CAR /CARG/' $pdb
sed -i .bu -e 's/LYS/LYP/' $pdb
sed -i .bu -e 's/HIS/HIE/' $pdb
sed -i .bu -e 's/CYS/CYN/' $pdb
sed -i .bu -e 's/CHI /CHIE/' $pdb

