#! /bin/bash

WEBB=webbgroup@pinotnoir.cm.utexas.edu
PRINTDIR=/Users/webbgroup/Desktop/Andrew/Printing

for arg in "$@" ; do

    if test "$arg" == "" ; then
        echo ; echo $'\a'Supply a file name to print, ie:
        echo print document.pdf ; echo
        exit
    fi

    new=`basename "$arg" | sed 's/ //g' | sed 's/^/temp./'`
    cp "$arg" $new

    if ! [[ $(file -b $new) == "PDF"* ]] ; then 
        echo ; echo "ALERT: $arg could not print"
        echo "-------> File must be in PDF format. "
        echo "-------> Please reprint document as 'print document.pdf'" ; echo ;  
## If the file type is not PDF, then print error message and do not Print
## If it is a PDF, proceed with printing procedure -- jtf 03 Mar 15
        exit
        fi
        
    scp $new $WEBB:$PRINTDIR
    echo; echo "Submitting $new to print..."
    ssh $WEBB 'lp -d Dell_3130cn_Color_Laser' $PRINTDIR/$new
    echo "Cleaning up..." ; echo
    ssh $WEBB rm $PRINTDIR/$new
    rm $new

done
