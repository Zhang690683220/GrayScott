module purge
module load gcc
module use /projects/community/modulefiles/
module load gcc/5.4/openmpi/3.1.2-kholodvl

# Start Dataspaces Server

cd ~/dataspaces/build_tcp/bin
rm -f conf
srun --mpi=pmix_v2 --mem-per-cpu=100000 --constraint=hal --ntasks=1 ./dataspaces_server -s 1 -c 32 >& log.server &

sleep 1s
while [ ! -f conf ]; do
    #echo "-- File conf is not yet available from server. Sleep more"
    sleep 1s
done
sleep 10s  # wait server to fill up the conf file

# Copy conf to working dir
cp ~/dataspaces/build_tcp/bin/conf ~/gray-scott/build_compression
cd ~/gray-scott/build_compression

# Run multiple tasks
srun --mpi=pmix_v2 --mem-per-cpu=20000 --constraint=slepner --ntasks=32 ./gray-scott ../simulation/settings-files.json