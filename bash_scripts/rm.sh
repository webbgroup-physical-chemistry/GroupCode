#! /bin/bash

if [ -z $1 ]; then
    echo $'\a'No file indicated for removal
    exit
    fi

until [ -z "$1" ]; do
    if [ -f $1 ]; then
        name=`basename $1`
        if [ -f $HOME/.Trash/$name ] ; then
            index=`ls $HOME/.Trash/$name* | wc -l`
            mv $1 $HOME/.Trash/$name.$index
        else
            mv $1 $HOME/.Trash/$name
            fi
        echo "$1 deleted"
    elif [ -d $1 ] ; then
        name=`basename $1`
        if [ -d $HOME/.Trash/$name ] ; then
            index=`ls $HOME/.Trash/$name* -d | wc -l`
            mv $1 $HOME/.Trash/$name.$index
        else
            mv $1 $HOME/.Trash/$name
            fi
        echo "$1 deleted"
    else
        echo $'\a'"File <$1> not found, nothing to delete..."
        fi
    shift
    done
