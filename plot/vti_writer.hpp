#ifndef __VTI_WRITER_H__
#define __VTI_WRITER_H__

#include <stdlib.h>
#include <stdio.h>
#include <type_traits>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>


#include <vtkImageDataGeometryFilter.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>

template <typename T>
void vti_writer(T *buffer, int nx, int ny, int nz, std::string filename)
{
    /* Write VTI file */
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(nx,ny,nz);
    if(std::is_same<T, int32_t>::value)
    {
        imageData->AllocateScalars(VTK_INT, 1);
    }
    else if(std::is_same<T, int64_t>::value)
    {
        imageData->AllocateScalars(VTK_LONG_LONG, 1);
    }
    else if(std::is_same<T, float>::value)
    {
        imageData->AllocateScalars(VTK_FLOAT, 1);
    }
    else if(std::is_same<T, double>::value)
    {
        imageData->AllocateScalars(VTK_DOUBLE, 1);
    }
    else
    {
        fprintf(stderr, "VTI_Writer doesn't support this data type\n");
        exit(1);
    }
    
    int* dims = imageData->GetDimensions();

    // Fill every entry of the image data 
    for (int z = 0; z < dims[2]; z++)
    {
        for (int y = 0; y < dims[1]; y++)
        {
            for (int x = 0; x < dims[0]; x++)
            {
                T* pixel = static_cast<T*>(imageData->GetScalarPointer(x,y,z));
                //pixel[0] = (float) x+512*y+512*512*z;
                pixel[0] = buffer[x+nx*y+nx*ny*z];
                //printf("%f ", pixel[0]);
            }
        }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    std::string filepath = "../img/";
    std::string suffix = ".vti";
    std::string path = filepath + filename +suffix;
    writer->SetFileName(path.c_str());
    writer->SetInputData(imageData);
    writer->Write();

    return;
}

//int vti_writer(int &buffer, int nx, int ny, int nz, std::string filename);
//int vti_writer(float &buffer, int nx, int ny, int nz, std::string filename);
//int vti_writer(double &buffer, int nx, int ny, int nz, std::string filename);



#endif