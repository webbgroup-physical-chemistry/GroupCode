#! /bin/bash

removed=$(find .. -name "*#" -print | xargs ls -l | awk '{ x += $5 } END {  x /= 1024 ; printf "Removed total %i KB\n", x }') 
find .. -name "*#" -print | xargs rm -v
echo $removed
