#!/bin/sh



# Init
module purge
module use /projects/community/modulefiles/
module load gcc/5.4/openmpi/3.1.2-kholodvl

# Start Dataspaces Server
cd ~/dataspaces/build_tcp/bin
rm -f conf
srun --mpi=pmix_v2 --mem-per-cpu=5000 --ntasks=1 ./dataspaces_server -s 1 -c 32 >& log.server &

wait