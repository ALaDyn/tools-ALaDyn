#BSUB -J aladyn
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -N
#BSUB -B
#BSUB -q hpc_inf
#BSUB -a openmpi
#BSUB -n 64
#BSUB -R "span[ptile=32]"
module load compilers/gcc-4.8.2
module load compilers/openmpi-1.8.1_gcc-4.8.2
/usr/share/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/mpirun.lsf env PSM_SHAREDCONTEXTS_MAX=8 /storage/gpfs_maestro/hpc/user/sinigardi/TNSA/NARA_EXPERIMENT/scan/ALaDyn >> opic.txt 2>> epic.txt
