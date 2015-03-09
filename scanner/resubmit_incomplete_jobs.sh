#! /bin/bash

for folder in ` find . -mindepth 1 -maxdepth 1 -type d `
do 
 cd $folder
 if [ ! -f spec04.dat ]
 then 
  echo "$folder da risottomettere"
#  qsub eurora-64.cmd
# else 
#  echo "$folder completata"
 fi
 cd ..
done

