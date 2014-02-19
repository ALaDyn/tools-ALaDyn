#! /bin/bash

if [ $# -lt 5 ]
then
 echo "In input devono esser passati:"
 echo "\$1 : NUMERO NODI BGQ RICHIESTI (nb: ciascuno contiene 16 CPU), min=64"
 echo "\$2 : NUMERO CPU PER NODO RICHIESTE"
 echo "\$3 : NUMERO ORE PER SINGOLA CPU RICHIESTE nel formato 'hh:mm:ss'"
 echo "\$4 : ACCOUNT DA ADDEBITARE"
 echo "\$5 : NOME JOB"
 echo "\$6 : opzionale: indirizzo email su cui inviare le notifiche"
 exit
fi

if [ $# == 6 ]
then
 EMAIL=$6
else
 USERNAME_SCRIPT=`whoami`
 HOST_SCRIPT="localhost"
 EMAIL="$USERNAME_SCRIPT@$HOST_SCRIPT"
fi


NUMBER_OF_CPUS=$(( $1 * $2 ))

 JOBFILE=job${NUMBER_OF_CPUS}.txt

 rm  ${JOBFILE}
 touch ${JOBFILE}
 chmod 755 ${JOBFILE}

 printf '#!/bin/bash\n' >> ${JOBFILE}
 printf '#\n' >> ${JOBFILE}
 printf '# @ account_no = %s\n' "$4" >> ${JOBFILE}
 printf '# @ job_name = %s\n' "$5" >> ${JOBFILE}
 printf '# @ output = opic.$(jobid)\n' >> ${JOBFILE}
 printf '# @ error = epic.$(jobid)\n' >> ${JOBFILE}
 printf '# @ shell = /bin/bash\n' >> ${JOBFILE}
 printf '# @ wall_clock_limit = %s\n' "$3" >> ${JOBFILE}
 printf '# @ job_type = bluegene\n' >> ${JOBFILE}
 printf '# @ notification = always\n' >> ${JOBFILE}
 printf '# @ bg_size = %s\n' "$1" >> ${JOBFILE}
 printf '# @ notify_user = %s\n' "${EMAIL}" >> ${JOBFILE}
 printf '# @ queue\n' >> ${JOBFILE}
 printf 'date\n' >> ${JOBFILE}
 printf 'runjob --np %s --ranks-per-node %s --exe ./rpic\n' "${NUMBER_OF_CPUS}" "$2" >> ${JOBFILE}
 printf 'date\n' >> ${JOBFILE}

