#! /bin/bash

#~ritchie/Utilities/python_scripts/NoSolIndex.py Rap_E30_K31+N27C.60-300.18.gro > nosol.ndx
#trjconv -f Rap_E30_K31+N27C.60-300.18.gro -o rap.gro -s Rap_E30_K31+N27C.60-300.5.tpr -n nosol.ndx
#trjconv -f Rap_E30_K31+N27C.60-300.18.xtc -o rap.xtc -s Rap_E30_K31+N27C.60-300.5.tpr -n nosol.ndx

# Make a tcl script!
# or use a single aligned template in which I know there are an equivalent number of CA

#echo 1 | pdb2gmx -f rap_aligned.pdb -o rap_aligned.gro -water none -p aligned -i aligned
#grompp -f /Volumes/Alfheim/ritchie/MD/x9_Umbrella/gromacs_mdp/diel_78_min.mdp -c rap_aligned.gro -p aligned -o rap_aligned.tpr
#echo C-Alpha System | trjconv -fit rot+trans -f rap.xtc -s rap_aligned.tpr 


