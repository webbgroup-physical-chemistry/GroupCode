#! /bin/bash

clear
for i in `seq 100`; do echo ; done
make clean
make
if [ -f g_CNC-HOH ] ; then
    clear
#   ./g_CNC-HOH -f test/rap_aligned.xtc -s test/rap_aligned.tpr -o test/hbonds.xvg -op test/persist.xvg -e 15 -oa test/residues.log
    time ./g_CNC-HOH -f test/production-1.xtc -s test/production-1.tpr -o test/hbonds.xvg -op test/persist.xvg -oa test/residues.log
#    ./g_CNC-HOH -f test/frame96.gro -s test/production-1.tpr -o test/blah.xvg


#    echo 0 1 | ./g_CNC-HOH -f test/rap_aligned.xtc -s test/rap_aligned.tpr -o test/blah.xvg
    rm test/*#
    fi
