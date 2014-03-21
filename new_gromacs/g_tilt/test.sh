#! /bin/bash

clear
for i in `seq 100`; do echo ; done
make clean
make
if [ -f g_tilt ] ; then
    clear
    echo 0 1 | ./g_tilt -f test/rap_aligned.gro -s test/rap_aligned.tpr -o test/blah.xvg -fr test/1LFD.tpr -n test/ral.ndx -v
#    echo 0 1 | ./g_tilt -f test/rap_aligned.xtc -s test/rap_aligned.tpr -o test/blah.xvg -fr test/1LFD.tpr -n test/ral.ndx
    rm test/*#
    fi
