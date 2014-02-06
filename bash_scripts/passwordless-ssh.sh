#! /bin/bash

target=$1
path=$2

if [ -z $target -o -z $path ] ; then 
	echo "usage: $0 <target> <path>"  
	echo "The <target> is computer@host"
	echo "The <path> is the location of target:/path/.ssh"
	exit
fi

if [ ! -f $HOME/.ssh/id_dsa.pub ]; then
	echo | ssh-keygen -t dsa
else
	echo "id_dsa.pug already exists, continuing..."
	fi
cp $HOME/.ssh/id_dsa.pub $HOME/.ssh/temp.key
scp $HOME/.ssh/temp.key $target:/$path/.ssh

ssh $target "cat /$path/.ssh/temp.key >> /$path/.ssh/authorized_keys"

ssh $target rm $path/.ssh/temp.key
rm $HOME/.ssh/temp.key

