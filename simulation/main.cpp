#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>


#include <mpi.h>

#include "../common/timer.hpp"
#include "gray-scott.h"

#include "dataspaces.h"
#include "zfp.h"


int ds_put(void* databuf, size_t bufsize, char* var_name, int timestep, const char* lock_name)
{
    // DataSpaces: Lock Mechanism
	// Usage: Prevent other process from modifying 
	// 	  data at the same time as ours
	// The NULL parameter is for a pointer to the 
	//MPI Communicator which we are not using in this example
	dspaces_lock_on_write(lock_name, NULL);

    sleep(2);

	//Name the Data that will be writen
	//char var_name[128];
	//sprintf(var_name, "ex_sample_data");

    printf("Timestep %d: put compressed data\n", timestep);
    

    // ndim: Dimensions for application data domain
	// In this case, our data array is 1 dimensional

	int ndim = 1;
    // Prepare LOWER and UPPER bound dimensions
	// In this example, we will put all data into a 
	// small box at the origin upper bound = lower bound = (0,0,0)
	// In further examples, we will expand this concept.
	uint64_t lb[3] = {0}, ub[3] = {0};

    ub[0] = bufsize;

    // DataSpaces: Put data array into the space
	// Usage: dspaces_put(Name of variable, version num, 
	// size (in bytes of each element), dimensions for bounding box,
	// lower bound coordinates, upper bound coordinates,
	// ptr to data buffer 
    //printf("debug1, bufsize = %d\n", bufsize);
	dspaces_put(var_name, timestep, 1*sizeof(char), ndim, lb, ub, databuf);
    dspaces_put_sync();
    // DataSpaces: Release our lock on the data
	dspaces_unlock_on_write(lock_name, NULL);

    return 0;

}

/*
void print_io_settings(const adios2::IO &io)
{
    std::cout << "Simulation writes data using engine type:              "
              << io.EngineType() << std::endl;
}
*/
void print_settings(const Settings &s)
{
    std::cout << "grid:             " << s.L << "x" << s.L << "x" << s.L
              << std::endl;
    std::cout << "steps:            " << s.steps << std::endl;
    std::cout << "plotgap:          " << s.plotgap << std::endl;
    std::cout << "F:                " << s.F << std::endl;
    std::cout << "k:                " << s.k << std::endl;
    std::cout << "dt:               " << s.dt << std::endl;
    std::cout << "Du:               " << s.Du << std::endl;
    std::cout << "Dv:               " << s.Dv << std::endl;
    std::cout << "noise:            " << s.noise << std::endl;
    std::cout << "output:           " << s.output << std::endl;
    std::cout << "adios_config:     " << s.adios_config << std::endl;
}

void print_simulator_settings(const GrayScott &s)
{
    std::cout << "process layout:   " << s.npx << "x" << s.npy << "x" << s.npz
              << std::endl;
    std::cout << "local grid size:  " << s.size_x << "x" << s.size_y << "x"
              << s.size_z << std::endl;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 1;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    MPI_Comm gcomm = MPI_COMM_WORLD;
    std::cout << "DEBUG1" << std::endl;
    dspaces_init(procs, 1, &gcomm, NULL);
    std::cout << "DEBUG2" << std::endl;

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Too few arguments" << std::endl;
            std::cerr << "Usage: gray-scott settings.json" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    Settings settings = Settings::from_json(argv[1]);

    GrayScott sim(settings, comm);
    sim.init();

    //adios2::ADIOS adios(settings.adios_config, comm, adios2::DebugON);
    //adios2::IO io_main = adios.DeclareIO("SimulationOutput");
    //adios2::IO io_ckpt = adios.DeclareIO("SimulationCheckpoint");

    //Writer writer_main(settings, sim, io_main);
    //Writer writer_ckpt(settings, sim, io_ckpt);

    //writer_main.open(settings.output);

    if (rank == 0) {
        //print_io_settings(io_main);
        std::cout << "========================================" << std::endl;
        print_settings(settings);
        print_simulator_settings(sim);
        std::cout << "========================================" << std::endl;
    }

#ifdef ENABLE_TIMERS
    Timer timer_total;
    Timer timer_compute;
    Timer timer_write;

#ifdef ENABLE_COMPRESSION
    Timer timer_compression;
#endif  

    std::ostringstream log_fname;
    log_fname << "gray_scott_pe_" << rank << ".log";
#ifdef ENABLE_COMPRESSION
    std::ofstream log(log_fname.str());
    log << "step\ttotal_gs\tcompute_gs\tcompression_gs\twrite_gs" << std::endl;
#else
    std::ofstream log(log_fname.str());
    log << "step\ttotal_gs\tcompute_gs\twrite_gs" << std::endl;
#endif    
#endif

    for (int i = 0; i < settings.steps;) {
#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_total.start();
        timer_compute.start();
#endif

        for (int j = 0; j < settings.plotgap; j++) {
            sim.iterate();
            i++;
        }

#ifdef ENABLE_TIMERS
        double time_compute = timer_compute.stop();
        MPI_Barrier(comm);
#endif

        if (rank == 0) {
            std::cout << "Simulation at step " << i
                      << " writing output step     " << i / settings.plotgap
                      << std::endl;
        }

        std::vector<double> u = sim.u_noghost();
        std::vector<double> v = sim.v_noghost();

#ifdef ENABLE_COMPRESSION

        double u_max =  *std::max_element(u.begin(), u.end());
        double u_min =  *std::min_element(u.begin(), u.end());
        double v_max =  *std::max_element(v.begin(), v.end());
        double v_min =  *std::min_element(v.begin(), v.end());

#endif

        double *u_array = (double*) malloc(u.size()*sizeof(double));
        double *v_array = (double*) malloc(v.size()*sizeof(double));

        std::copy(u.begin(), u.end(), u_array);
        std::copy(v.begin(), v.end(), v_array);

#ifdef ENABLE_COMPRESSION

#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_compression.start();
#endif



        /* Compression Init */
        zfp_type type;     /* array scalar type */
        zfp_field *u_field, *v_field;  /* array meta data */
        zfp_stream *u_zfp, *v_zfp;   /* compressed stream */
        void *u_buffer, *v_buffer;      /* storage for compressed stream */
        size_t u_bufsize, v_bufsize;    /* byte size of compressed buffer */
        bitstream *u_stream, *v_stream; /* bit stream to write to or read from */
        size_t u_zfpsize, v_zfpsize;    /* byte size of compressed stream */
        double rate = 32;   /* bit size of each compressed double */
        double torlerance = 1e-3;

        /* allocate meta data for the 3D array a[nz][ny][nx] */
        type = zfp_type_double;
        u_field = zfp_field_3d(u_array, type, sim.size_x, sim.size_y, sim.size_z);
        v_field = zfp_field_3d(v_array, type, sim.size_x, sim.size_y, sim.size_z);

        /* allocate meta data for a compressed stream */
        u_zfp = zfp_stream_open(NULL);
        v_zfp = zfp_stream_open(NULL);

        /* set compression mode and parameters via one of three functions */
        /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */ 
        /*  zfp_stream_set_precision(zfp, precision); */
        zfp_stream_set_accuracy(u_zfp, (u_max-u_min)*torlerance);
        zfp_stream_set_accuracy(v_zfp, (v_max-v_min)*torlerance);

        /* allocate buffer for compressed data */
        u_bufsize = zfp_stream_maximum_size(u_zfp, u_field);
        v_bufsize = zfp_stream_maximum_size(v_zfp, v_field);
        u_buffer = malloc(u_bufsize);
        v_buffer = malloc(v_bufsize);

        /* associate bit stream with allocated buffer */
        u_stream = stream_open(u_buffer, u_bufsize);
        v_stream = stream_open(v_buffer, v_bufsize);
        zfp_stream_set_bit_stream(u_zfp, u_stream);
        zfp_stream_set_bit_stream(v_zfp, v_stream);
        zfp_stream_rewind(u_zfp);
        zfp_stream_rewind(v_zfp);

        /* compress array and output compressed stream */
        u_zfpsize = zfp_compress(u_zfp, u_field);
        v_zfpsize = zfp_compress(v_zfp, v_field);

#ifdef ENABLE_TIMERS
        double time_compression = timer_compression.stop();
        MPI_Barrier(comm);
#endif        
        
        if (!(u_zfpsize && v_zfpsize)) {
            fprintf(stderr, "compression failed\n");
            exit(1);
        }
        else {
#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_write.start();
#endif


            char u_var_name[128];
            char v_var_name[128];
		    sprintf(u_var_name, "u_compressed_data");
            sprintf(v_var_name, "v_compressed_data");
            const char* u_lock_name = "u_compressed_lock";
            const char* v_lock_name = "v_compressed_lock";
            ds_put(u_buffer, u_zfpsize, u_var_name, i, u_lock_name);
            ds_put(v_buffer, v_zfpsize, v_var_name, i, v_lock_name);

        }

        /* clean up */
        zfp_field_free(u_field);
        zfp_field_free(v_field);
        zfp_stream_close(u_zfp);
        zfp_stream_close(v_zfp);
        stream_close(u_stream);
        stream_close(v_stream);
        free(u_buffer);
        free(v_buffer);
#else

#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_write.start();
#endif
        char u_var_name[128];
        char v_var_name[128];
		sprintf(u_var_name, "u_data");
        sprintf(v_var_name, "v_data");
        const char* u_lock_name = "u_lock";
        const char* v_lock_name = "v_lock";
        ds_put(u_array, u.size()*sizeof(double), u_var_name, i, u_lock_name);
        ds_put(v_array, v.size()*sizeof(double), v_var_name, i, v_lock_name);


#endif

        //writer_main.write(i, sim);
        /*
        if (settings.checkpoint &&
            i % (settings.plotgap * settings.checkpoint_freq) == 0) {
            writer_ckpt.open(settings.checkpoint_output);
            writer_ckpt.write(i, sim);
            writer_ckpt.close();
        }
        */
#ifdef ENABLE_TIMERS
        double time_write = timer_write.stop();
        double time_step = timer_total.stop();
        MPI_Barrier(comm);

#ifdef ENABLE_COMPRESSION
        log << i << "\t" << time_step << "\t" << time_compute << "\t"
            << time_compression << "\t" << time_write << std::endl;
#else
        log << i << "\t" << time_step << "\t" << time_compute << "\t"
            << time_write << std::endl;
#endif
#endif

    free(u_array);
    free(v_array);
    }

    //writer_main.close();

#ifdef ENABLE_TIMERS

#ifdef ENABLE_COMPRESSION
    log << "total\t" << timer_total.elapsed() << "\t" << timer_compute.elapsed()
        << "\t" << timer_compression.elapsed() << "\t" << timer_write.elapsed() << std::endl;
#else
    log << "total\t" << timer_total.elapsed() << "\t" << timer_compute.elapsed()
        << "\t" << timer_write.elapsed() << std::endl;
#endif
    log.close();
#endif

    MPI_Finalize();
}
