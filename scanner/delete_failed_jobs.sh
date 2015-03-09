#! /bin/bash

for folder in ` find . -mindepth 1 -maxdepth 1 -type d `
do 
 cd $folder
 if [ -s epic.txt ]
 then 
  echo "$folder fallita"
  cd ..
  rm -rf $folder
 else
  cd ..
 fi
done

