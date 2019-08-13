


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkPoints.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>

#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkCurvatures.h>
#include <vtkPointData.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>

#include <vtkDoubleArray.h>

#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>

#include <vtkMetaImageWriter.h>
#include <vtkPolyDataWriter.h>

#include <vtkImageDilateErode3D.h>

#include <vnl/vnl_vector.h>

#include <map>
#include <vector>

#define PI 3.14159265

using namespace std;

void help(char* execname){
    cout<<"Computes a binary image from the wm surface and the gm surface"<<endl;
    cout<<"Usage: "<<execname<<" -d <patients dir> -p <patient name> [options]"<<endl;
    cout<<"options:"<<endl;
    cout<<"--h --help show help menu"<<endl;    

}


void PolyDataImageStencilExport(vtkPolyData* polydata, vtkImageData *imagein, double* bounds){



    //vtkPolyData* polydata = (vtkPolyData*) poly->GetInput();

    vtkSmartPointer<vtkPolyDataToImageStencil> polytostencil = vtkSmartPointer<vtkPolyDataToImageStencil>::New();

    polytostencil->SetInput(polydata);

    polytostencil->Update();



    //double *bounds = polydata->GetBounds();

   // vtkSmartPointer<vtkImageData> imagein = vtkSmartPointer<vtkImageData>::New();



    imagein->SetExtent(bounds[0] - 1, bounds[1] + 1, bounds[2] - 1, bounds[3] + 1, bounds[4] - 1, bounds[5] + 1);

    imagein->SetScalarTypeToUnsignedShort();

    imagein->AllocateScalars();





    int* extent = imagein->GetExtent();



    for (int x = extent[0]; x <= extent[1]; x++){

        for (int y = extent[2]; y <= extent[3]; y++){

            for (int z  =extent[4]; z <= extent[5]; z++){

                unsigned short* pixel = static_cast<unsigned short*>(imagein->GetScalarPointer(x,y,z));

                *pixel = 0;

            }

        }

    }



    vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();

    stencil->SetInput(imagein);

    stencil->SetStencil(polytostencil->GetOutput());

    stencil->ReverseStencilOn();

    stencil->SetBackgroundValue(128);

    stencil->Update();

    imagein->DeepCopy(stencil->GetOutput());



}


int main(int argc, char *argv[])
{


    if (argc == 1){
        help(argv[0]);
        return 0;
    }

    string dirname = "";
    vector< std::string > patientname;

    for(int i = 1; i < argc; i++){

        if(string(argv[i]) == "--h" || string(argv[i]) == "--help"){
            help(argv[0]);
            return 0;
        }else if (string(argv[i]) == "-d"){
            if(i + 1 >= argc){
                std::cout << "Patient directory missing"<< std::endl ;
                return EXIT_FAILURE;
            }
            dirname = argv[i + 1];


        }else if(string(argv[i]) == "-p"){

            //imagesize[0] = atoi(argv[i+1]);
            patientname.push_back(argv[i + 1]);

        }

    }

    if(patientname.size() == 0){
        cout<<"Error Patient name missing!"<<endl;
        return EXIT_FAILURE;
    }

    cout<<"patient: "<<patientname[0]<<endl;


    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();

    string origfilename = dirname;
    origfilename.append("/");
    origfilename.append(patientname[0]);
    origfilename.append("/surf/vtk/lh.orig.vtk");

    reader->SetFileName(origfilename.c_str());
    reader->Update();
    vtkPolyData* orig = reader->GetOutput();

    string pialfilename = dirname;
    pialfilename.append("/");
    pialfilename.append(patientname[0]);
    pialfilename.append("/surf/vtk/lh.pial.vtk");

    vtkSmartPointer<vtkPolyDataReader> reader2 =  vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(pialfilename.c_str());
    reader2->Update();
    vtkPolyData* pial = reader2->GetOutput();
    double *bounds = pial->GetBounds();


    vtkSmartPointer<vtkImageData> pialimage = vtkSmartPointer<vtkImageData>::New();
    PolyDataImageStencilExport(pial, pialimage, bounds);


    vtkSmartPointer<vtkImageData> origimage = vtkSmartPointer<vtkImageData>::New();
    PolyDataImageStencilExport(orig, origimage, bounds);


    int* extent = pialimage->GetExtent();

    for (int z  =extent[4] + 1; z < extent[5]; z++){
        for (int x = extent[0] + 1; x < extent[1]; x++){
            for (int y = extent[2] + 1; y < extent[3]; y++){

                unsigned short* pialpixel = static_cast<unsigned short*>(pialimage->GetScalarPointer(x,y,z));
                unsigned short* origpixel = static_cast<unsigned short*>(origimage->GetScalarPointer(x,y,z));


                if(*pialpixel != 0 && *origpixel != 0){
                    *pialpixel = 0;
                }
            }
        }
    }


    vtkSmartPointer<vtkMetaImageWriter> imagewriter = vtkSmartPointer<vtkMetaImageWriter>::New();
    imagewriter->SetInput(pialimage);

    string imagefile = dirname;
    imagefile.append("/");
    imagefile.append(patientname[0]);
    imagefile.append("/surf/lh.binimage.mhd");
    imagewriter->SetFileName(imagefile.c_str());
    imagewriter->Write();

    return 0;
}


