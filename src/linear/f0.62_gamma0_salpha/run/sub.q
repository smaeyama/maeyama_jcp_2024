#!/bin/bash

###  NOTE  ###
###  Plasma simulator, NEC SX-Aurora TSUBASA A412-8 (NIFS, 2020)
###
###  - Computation nodes (total 4320 VE (Vector engine))
###      VE model: Type 10AE (8cores)
###      Peak performance: DP 2.433 TFLOPS per VE
###      Memory: HBM2 48 GiB
###      Memory Bandwidth: ? GB/s per node
###
###      (For now, flat MPI is recommended.)
###
###  - Interconnect
###      Infiniband HDR200 x2, 1000BASE-Tx1, BMC
###
###  - Job class : Computation server (SX-Aurora)
###      small     :    1 - 16   VE, 15 min., 1 run/ 1 submit
###      small24VE :    1 - 4    VE, 24 hour, 8 run/16 submit
###      small24VH :    8 - 32   VE, 24 hour, 8 run/16 submit
###      medium    :   40 - 768  VE, 10 hour, 4 run/ 8 submit
###      large     : 1920 - 2160 VE, 10 hour, 1 run/ 4 submit
###      large1h   : 1920 - 2160 VE,  1 hour, 1 run/ 2 submit
###      debug     :    8 - 16   VE, 30 min., 1 run/ 1 submit, interactive
###
###  - Job class : Data analysis server (LX)
###      gpu-b : 1 - 4 Servers, 10 hour, 1 run/2 submit
###      gpu-i : 1 - 2 Servers, 10 hour, 1 run/1 submit, interactive
###
###  - Commands
###      (Submit a batch job : "qsub sub.q") Use shoot script for GKV.
###      Check job status    : "qstat -a"
###      Delete job          : "qdel JOBID"
###      Show budget info    : "pstime"
###      Show disk usage     : "lsquota"
##############

#PBS -q small24VE                 # queue name
#PBS --group=22265                # resource group
#PBS -T necmpi                    # necessary for MPI job
#PBS -l elapstim_req=00:15:00     # elapsed time limit

#PBS --venode=2                   # total number of VE
#PBS --venum-lhost=2              # number of VE per a logical node
#### --venum-lhost=8              # number of VE per a logical node
#PBS -v OMP_NUM_THREADS=1         # number of threads per MPI process

MPI_procs=16                      # number of MPI processes (= venode*8 for flat MPI)

#PBS -v VE_FORT_SETBUF=10240
#PBS -v FTRACE=YES
#PBS -v NMPI_PROGINF=DETAIL
#PBS -v NMPI_SEPSELECT=3

#PBS -v LANG=C

source /etc/profile.d/modules.sh

module load NECNLC-sx
# module load NECNLC-mpi-sx
### For NetCDF
module load netcdf-parallelIO-fortran-sx


### Working directory 
DIR=%%DIR%%
LDM=gkvp.exe
NL=gkvp_namelist.%%%


date
cd ${DIR}
export fu05=${DIR}/${NL}


#cat << 'EOF-S' > ./mpisep.sh
##!/bin/sh
#ulimit -s unlimited
#ID=${MPIUNIVERSE}.`printf "%05d" ${MPIRANK}`
#case ${NMPI_SEPSELECT:-${MPISEPSELECT:-2}} in
#1) exec $* 1>> stdout.${ID}                  ;;
#2) exec $*                  2>> stderr.${ID} ;;
#3) exec $* 1>> stdout.${ID} 2>> stderr.${ID} ;;
#4) exec $* 1>> std.${ID}    2>&1             ;;
#*) exec $* ;;
#esac
#EOF-S
#chmod 777 ./mpisep.sh
#
##---( time mpiexec -v -nn ${_NECMPI_VH_NUM_NODES} -ve 0-7 -ppn 64 ./mpisep.sh ./${LDM} ) > log.mpi 2>&1
#( time mpiexec -v -nn ${_NECMPI_VH_NUM_NODES} -ve 0-7 -ppn 64 -n ${MPI_procs} ./mpisep.sh ./${LDM} ) > log.mpi 2>&1

mpirun -n ${MPI_procs} ${DIR}/${LDM}

date
#touch complete


#---#PBS -l coresz_prc=10
#---#PBS --venum-lhost=8
#---#PBS -b 4                     # number of nodes

