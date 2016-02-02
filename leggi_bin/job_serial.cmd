#!/bin/bash

# @ account_no = my_cineca_computing_account
# @ shell = /bin/bash
# @ job_type = serial
# @ job_name = postproc.$(jobid)
# @ output = postproc.out.$(jobid)
# @ error = postproc.err.$(jobid)
# @ wall_clock_limit = 0:30:00
# @ class = serial
# @ notification = always
# @ notify_user = myemail@address
# @ queue

/path/to/my/executable argumentlist -goes here

