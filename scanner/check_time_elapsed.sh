#!/bin/bash

declare -a FOLDERS
FOLDERS=`find . -mindepth 1 -maxdepth 1 -type d`

tempo_max=0
for folder in $FOLDERS
do
 cd $folder
  tempo=`tail -2 opic.txt |head -1`
  tempo=${tempo#*  }
  tempo=${tempo%,*}
  tempo_int=${tempo%.*}
  if [[ $tempo_int -gt $tempo_max ]]
  then 
   tempo_max=$tempo_int
  fi
  echo $folder : $tempo
 cd ..
done
echo
tempo_max_h=`echo " $tempo_max / 3600 " | bc -l`
printf "Slowest required %d seconds, equivalent to %3.1f hours\n" "$tempo_max" "$tempo_max_h"


