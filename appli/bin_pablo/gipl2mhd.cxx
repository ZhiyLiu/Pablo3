//----------------------------------------------------------------------
//	File:			gipl2mhd.cxx
//  convert a giplimage to mhd
//  Author:       Prieto
//----------------------------------------------------------------------





#include "itkGiplImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


using namespace std;					// make std:: accessible

void help(char* execname){
    cout<<"convert a giplimage to mhd"<<endl;
    cout<<"Usage: "<<execname<<" <giplimage>"<<endl;
    cout<<"--h --help show help menu"<<endl;

}


int main(int argc, char **argv)
{


    if (argc < 2)
    {
        help(argv[0]);
        return EXIT_FAILURE;
    }


    std::string inputFilename = argv[1];


    typedef unsigned short PixelType;
    typedef itk::Image< PixelType, 3 > InputImageType;

    typedef itk::ImageFileReader< InputImageType > ImageFileReaderType;
    ImageFileReaderType::Pointer imagereader = ImageFileReaderType::New();

    imagereader->SetFileName(inputFilename.c_str());
    imagereader->Update() ;



    string filenamefeature = inputFilename.substr(0, inputFilename.find_last_of("."));
    filenamefeature.append(".mhd");


    typedef itk::ImageFileWriter< InputImageType > ImageWriterType;
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName(filenamefeature.c_str());
    writer->SetInput(  imagereader->GetOutput() );
    writer->Update();





    return 0;
}



