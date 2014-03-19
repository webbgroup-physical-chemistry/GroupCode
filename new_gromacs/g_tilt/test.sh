#! /bin/bash

clear
make clean
make
if [ -f g_tilt ] ; then
    clear
#    ./g_tilt -f test/Rap_E30_K31+N27C.60-300.18.gro -s test/Rap_E30_K31+N27C.60-300.5.tpr -o test/blah.xvg -fr test/1LFD.tpr -n index.ndx
    ./g_tilt -f test/rap_aligned.gro -s test/rap_aligned.tpr -o test/blah.xvg -fr test/1LFD.tpr -n test/ral.ndx
    fi  
