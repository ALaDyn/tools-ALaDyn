#!/bin/bash
# @ account_no = Pra05_0994
# @ shell  = /bin/bash
# @ job_type = serial
# @ job_name =postproc.$(jobid)
# @ output = postproc.out.$(jobid)
# @ error = postproc.err.$(jobid)
# @ wall_clock_limit = 0:30:00
# @ class = serial
# @ notification = always
# @ notify_user = stefano.sinigardi@gmail.com
# @ queue

date
/fermi/home/userexternal/ssinigar/bin/leggi_binario Prpout10 -noswap 
date
