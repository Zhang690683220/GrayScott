#include "vti_writer.h"

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
                int* pixel = static_cast<int*>(imageData->GetScalarPointer(x,y,z));
                //pixel[0] = (float) x+512*y+512*512*z;
                pixel[0] = buffer[x+512*y+512*512*z];
                //printf("%f ", pixel[0]);
            }
        }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(imageData);
    writer->Write();

    return 0;
}



