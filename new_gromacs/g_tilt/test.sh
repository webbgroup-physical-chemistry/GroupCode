#! /bin/bash

clear
for i in `seq 100`; do echo ; done
#make clean
make
if [ -f g_tilt ] ; then
    clear
    echo 0 1 | ./g_tilt \
-f test_files/rap_aligned.gro \
-s test_files/rap_aligned.tpr \
-o test_files/blah.xvg \
-fr test_files/1LFD.tpr \
-n test_files/ral.ndx \
-v
time \
echo 0 1 | ./g_tilt \
-f test_files/rap_aligned.xtc \
-s test_files/rap_aligned.tpr \
-o test_files/blah_xtc.xvg \
-fr test_files/1LFD.tpr \
-n test_files/ral.ndx \

    rm test_files/*#
    fi
