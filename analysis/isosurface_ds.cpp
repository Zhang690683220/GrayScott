#include <iostream>
#include <sstream>

#include <vtkAppendPolyData.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>



#include "../common/timer.hpp"
#include "../simulation/settings.h"
#include "analysis.h"
#include "reader.h"

vtkSmartPointer<vtkPolyData>
compute_isosurface(const Analysis &anly, const std::vector<double> &field, double isovalue)
{
    // Convert field values to vtkImageData
    auto importer = vtkSmartPointer<vtkImageImport>::New();
    importer->SetDataSpacing(1, 1, 1);
    importer->SetDataOrigin(static_cast<double>(anly.offset_x), static_cast<double>(anly.offset_y),
                            static_cast<double>(anly.offset_z));
    importer->SetWholeExtent(0, static_cast<int>(anly.size_x) - 1, 0,
                             static_cast<int>(anly.size_y) - 1, 0,
                             static_cast<int>(anly.size_z) - 1);
    importer->SetDataExtentToWholeExtent();
    importer->SetDataScalarTypeToDouble();
    importer->SetNumberOfScalarComponents(1);
    importer->SetImportVoidPointer(const_cast<double *>(field.data()));

    // Run the marching cubes algorithm
    auto mcubes = vtkSmartPointer<vtkMarchingCubes>::New();
    mcubes->SetInputConnection(importer->GetOutputPort());
    mcubes->ComputeNormalsOn();
    mcubes->SetValue(0, isovalue);
    mcubes->Update();

    // Return the isosurface as vtkPolyData
    return mcubes->GetOutput();
}

void write_vtk(const std::string &fname,
               const vtkSmartPointer<vtkPolyData> polyData,
               int rank, int timestep)
{
    std::string filename = fname + ".Rank." + std::to_string(rank) + ".Time." + std::to_string(timestep) + ".vtk";
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fname.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}


void print_settings(const Settings &s)
{
    std::cout << "grid:             " << s.L << "x" << s.L << "x" << s.L
              << std::endl;
    std::cout << "steps:            " << s.steps << std::endl;
    std::cout << "plotgap:          " << s.plotgap << std::endl;
    std::cout << "output:           " << s.output << std::endl;
}

void print_simulator_settings(const Analysis &a)
{
    std::cout << "process layout:   " << a.npx << "x" << a.npy << "x" << a.npz
              << std::endl;
    std::cout << "local grid size:  " << a.size_x << "x" << a.size_y << "x"
              << a.size_z << std::endl;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 5;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    if (argc < 3) {
        if (rank == 0) {
            std::cerr << "Too few arguments" << std::endl;
            std::cerr << "Usage: isosurface_ds settings.json isovalues..." << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    Settings settings = Settings::from_json(argv[1]);

    std::vector<double> isovalues;
    for (int i = 3; i < argc; i++) {
        isovalues.push_back(std::stod(argv[i]));
    }

    Analysis anly(settings, comm);
    anly.init();

    if (rank == 0) {
        std::cout << "========================================" << std::endl;
        print_settings(settings);
        print_simulator_settings(anly);
        std::cout << "========================================" << std::endl;
    }

    Reader dsreader(comm, procs, 2);

    for (int i = 0; i < settings.steps;) {

        i+=settings.plotgap;

        if (rank == 0) {
            std::cout << " reading input step     " << i
                      << "Isosurface at step " << i / settings.plotgap 
                      << std::endl;
        }

        dsreader.read(anly, i);

        auto appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

        for (const auto isovalue : isovalues) {
            auto polyData = compute_isosurface(anly, anly.u_noghost(), isovalue);
            appendFilter->AddInputData(polyData);
        }

        appendFilter->Update();

        std::string fname = "grayscott_isosurface";
        write_vtk(fname, appendFilter->GetOutput(),rank, i);

        
    }

    MPI_Finalize();
}