#ifndef __READER_H__
#define __READER_H__

#include <mpi.h>
#include <string>
#include <vector>

#include "../simulation/settings.h"
#include "analysis.h"
#include "dataspaces.h"

class Reader
{
public:
    Reader(MPI_Comm comm ,int proc, int appid);
    
    void read(const Analysis &analysis, int step);
    void close() { dspaces_finalize();}

private:
    MPI_Comm &comm;
};

#endif
