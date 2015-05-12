#! /bin/bash
module load profile/advanced
module load gnu/4.9.2
module load boost/1.57.0--gnu--4.9.2
g++ -I/cineca/prod/libraries/boost/1.57.0/gnu--4.9.2/include/ -O2 -std=c++11 -o leggi_diag leggi_diag.cpp

