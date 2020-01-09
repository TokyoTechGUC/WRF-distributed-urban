#!/bin/bash
#$ -cwd
#$ -l f_node=6
#$ -l h_rt=01:30:00
#$ -N present_hcmc_08_landscan_urb_param_correct_ahe
#$ -m abe
#$ -M do.k.aa@m.titech.ac.jp

. /etc/profile.d/modules.sh
. $SGE_O_HOME/.bashrc

mpirun -ppn 26 ./wrf.exe
