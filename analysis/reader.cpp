#include "reader.h"
#include <iostream>



Reader::Reader(MPI_Comm comm ,int proc, int appid)
    :comm(comm)
{

     dspaces_init(proc, appid, &comm, NULL);
}

void Reader::read(const Analysis &analysis, int step)
{


    std::vector<double> u = analysis.u_noghost();
    std::vector<double> v = analysis.v_noghost();

    std::string VarNameU = "grayscott_u";
    std::string VarNameV = "grayscott_v";

    uint64_t lb[3] = {0}, ub[3] = {0};

    lb[0] = sim.offset_x;
    lb[1] = sim.offset_y;
    lb[2] = sim.offset_z;

    ub[0] = sim.offset_x+sim.size_x-1;
    ub[1] = sim.offset_y+sim.size_y-1;
    ub[2] = sim.offset_z+sim.size_z-1;

    MPI_Barrier(comm);

    std::string LockNameU = VarNameU + "_lock";
    std::string LockNameV = VarNameV + "_lock";

    dspaces_lock_on_read(LockNameU.c_str(), &comm);
    int status = dspaces_get(VarNameU.c_str(),step, sizeof(double) ,3 ,lb ,ub ,u.data());
    dspaces_unlock_on_read(LockNameU.c_str(), &comm);

    dspaces_lock_on_read(LockNameV.c_str(), &comm);
    status = dspaces_get(VarNameV.c_str(),step, sizeof(double) ,3 ,lb ,ub ,v.data());
    dspaces_unlock_on_read(LockNameV.c_str(), &comm);

    if (status!=0){
        std::cout << "error for the ts get at step " << step << "with status " << status << std::endl;
    }
}