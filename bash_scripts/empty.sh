#! /bin/bash

ls $HOME/.Trash/
spaceused=`du -sh $HOME/.Trash | awk '{print $1}'`
echo ; echo "Currently $spaceused in .Trash"
echo ; echo "WARNING! YOU ARE ABOUT TO DELETE EVERYTHING IN THE TRASH FILE," ; echo "ARE YOU SURE ABOUT THIS? THERE IS NO GOING BACK... " ; echo
read answer
if [ ${answer:0:1} = "y" -o ${answer:0:1} = "Y" ] ; then
    echo ; echo "ARE YOU ABSOLUTELY POSITIVE?  THIS IS YOUR LAST WARNING..." ; echo
    read answer
    if [ ${answer:0:1} = "y" -o ${answer:0:1} = "Y" ] ; then
        echo ; echo Deleting...
       rm -rf $HOME/.Trash/* # $HOME/.Trash/.*
        echo ; echo Trash has been emptied ; echo 
        df -h | awk '/home/{print}'
    else 
        echo ; echo Trash has not been emptied ; echo
        fi
else
    echo ; echo Trash has not been emptied ; echo 
    fi
