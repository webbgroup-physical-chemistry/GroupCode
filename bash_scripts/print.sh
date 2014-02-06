#! /bin/bash


for arg in "$@" ; do

    if test "$arg" == "" ; then
        echo ; echo $'\a'Supply a file name to print, ie:
        echo print document.pdf ; echo
        exit
    fi

    new=`basename "$arg" | sed 's/ //g' | sed 's/^/temp./'`
    cp "$arg" $new
    scp $new webbgroup@pinotnoir.cm.utexas.edu:/Users/webbgroup/Desktop/Andrew/Printing
    echo; echo "Submitting $new to print..."
    ssh webbgroup@pinotnoir.cm.utexas.edu 'lp -d Dell_3130cn_Color_Laser' /Users/webbgroup/Desktop/Andrew/Printing/$new
    echo "Cleaning up..." ; echo
    ssh webbgroup@pinotnoir.cm.utexas.edu rm /Users/webbgroup/Desktop/Andrew/Printing/$new
    rm $new

done