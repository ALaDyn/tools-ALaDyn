#! /bin/bash

JOB_LIST="my_jobs.txt"

cat ${JOB_LIST} | while read line
do 
JOB_NAME=$(echo "$line" | awk '{print $1}')
qdel "${JOB_NAME}"
done

