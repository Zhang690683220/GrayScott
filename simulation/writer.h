#ifndef __WRITER_H__
#define __WRITER_H__

#include <mpi.h>
#include <string>
#include <vector>

#include "gray-scott.h"
#include "settings.h"
#include "dataspaces.h"

class Writer
{
public:
    Writer(MPI_Comm comm ,int proc, int appid);
    
    void write(int step, const GrayScott &sim);
    void write(const GrayScott &sim, MPI_Comm comm, int step);
    void close() { dspaces_finalize();}

protected:
    std::string m_serverMasterAddr;
};

#endif
