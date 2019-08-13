
#define PI 3.14159265

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
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkQuad.h>


#include <vtkPointLocator.h>
#include <vtkSphereSource.h>
#include <vtkbestcircle.h>
#include <vtkAxesActor.h>
#include <vtkMetaImageReader.h>

#include <vtksrep.h>
#include <vtksrepinterpolatemedialsheet.h>
#include <vtksrepvisuprimitives.h>
#include <vtksrepinterpolatemedialcrestcurve.h>
#include <vtksrepinterpolatemedialspokes.h>
#include <vtksrepinterpolatecrestspokes.h>
#include <vtksrepinterpolatecrestspokesquartic.h>

#include <vtkoptimizeCortexAtomsDis.h>
#include <vtkoptimizeCortexAtomsCurv.h>

#include <vtkoptimizeSRepSampling.h>

#include "vnl/algo/vnl_levenberg_marquardt.h"
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtksrepinterpolatemedialspokeshermite.h>
#include <vnl/vnl_vector.h>
#include "vtkinterpolatecurve.h"


#include "vtkBYUWriter.h"
#include "vtkquadmeshtotriangularmesh.h"
#include "vtksreptogridlayout.h"

#include <map>


#include "P3DControl.h"
#include "ControlParms.h"

using namespace std;


double m_Min = 0.01;
double m_Max = 0.5;
int m_NumP = 40;
int m_Inter = 3;
bool m_RenderSphere = 0;
bool m_StartInteractor = 1;
bool m_Optimize = 1;
string m_Hemi = "lh";
int m_Surf = 0;

ControlParms * globalControl;	// Read the user's preferences file
int globalVerbosity;			// Current verbosity level of Pablo
P3DControl* p3d;

M3DQuadFigure* NewQuadFigure(const int numr, int numc){
    globalControl = new ControlParms(NULL, -1, false);	// Ignore user's preferences
    globalVerbosity = 0;
    globalControl->setDefault(OutputVerbosity, globalVerbosity);
    globalControl->setDefault(ReorderModels, 0);
    globalControl->setDefault(SmoothImages, false);
    globalControl->setDefault(ConvertImages, false);
    globalControl->setDefault(ByteOrder, 1);
    globalControl->setDefault(CompressImages, true);
    globalControl->setDefault(ShowLandmarks, true);

    p3d = new P3DControl(10);
    float color[3] = {0.2,0.5,0.8};
    p3d->addQuadFigure(numr, numc, "Cortex", color);

    M3DObject* m3dobject = p3d->getObjectPtr();
    M3DFigure * figure = m3dobject->getFigurePtr( 0 ) ;
    //figure->setLandmarkName(0,0)
   return dynamic_cast<M3DQuadFigure*>( figure );
}

void SaveQuadFigure(M3DQuadFigure* fig, const char* filename){
    p3d->write(filename);
}

void help(char* execname){
    cout<<"Registration using the structure with least points in the group of subjects"<<endl;
    cout<<"The stereographic projection is done for each subject and an image of the unfolded cortex"<<endl;
    cout<<"is created. This image is use to produce the registration using elastix."<<endl;
    cout<<"allLabelsRegistration.txt is created, the file containing the moving image for all the cases."<<endl;
    cout<<"Usage: "<<execname<<" -d <patients dir> -pfn <patient file name> -elastixp <path to elastix> [options]"<<endl;
    cout<<"options:"<<endl;
    cout<<"--h --help show help menu"<<endl;
    cout<<"-pfn <patient filename> list of file with name of patients <one per line>"<<endl;
    cout<<"-elastixp <elastixp> path to elastix"<<endl;


}

void createGridPoints(vtkSmartPointer< vtkImageData > img, int label, int count, bool all, string filename){
    int extent[6];
    img->GetExtent(extent);
    double origin[3];
    img->GetOrigin(origin);

    FILE * fgrid;
    fgrid = fopen (filename.c_str(),"w");
    if (fgrid!=NULL)
    {
        fputs ("point\n",fgrid);
        char buff[50];        
        sprintf(buff, "%d\n", count);
        fputs (buff,fgrid);

        for(int i = extent[0]; i <= extent[1]; i++){
            for(int j = extent[2]; j <= extent[3]; j++){
                for(int k = extent[4]; k <= extent[5]; k++){
                    float* ptr = (float*)img->GetScalarPointer(i, j, k);
                    if(*ptr == label && !all){
                        float x = (float)i + origin[0];
                        float y = (float)j + origin[1];
                        sprintf(buff, "%f %f\n", x, y);
                        fputs (buff,fgrid);
                    }else if(all){
                        float x = (float)i + origin[0];
                        float y = (float)j + origin[1];
                        sprintf(buff, "%f %f\n", x, y);
                        fputs (buff,fgrid);
                    }
                }
            }
        }

        fclose (fgrid);
    }



}

void createMask(vtkSmartPointer< vtkImageData > img, int label, string filename){

    vtkSmartPointer< vtkImageData > outmask = vtkSmartPointer< vtkImageData >::New();
    outmask->DeepCopy(img);

    int extent[6];
    outmask->GetExtent(extent);

    for(int i = extent[0]; i <= extent[1]; i++){
        for(int j = extent[2]; j <= extent[3]; j++){
            for(int k = extent[4]; k <= extent[5]; k++){
                float* ptr = (float*)outmask->GetScalarPointer(i, j, k);
                if(*ptr != label){
                    *ptr = 0;
                }else{
                    *ptr = 1;
                }
            }
        }
    }

    vtkSmartPointer< vtkMetaImageWriter > writer = vtkSmartPointer< vtkMetaImageWriter >::New();
    writer->SetInput(outmask);
    writer->SetFileName(filename.c_str());
    writer->Write();
}

int countLabels(vtkSmartPointer< vtkImageData > img, int label){

    int extent[6];
    img->GetExtent(extent);
    int count = 0;
    for(int i = extent[0]; i <= extent[1]; i++){
        for(int j = extent[2]; j <= extent[3]; j++){
            for(int k = extent[4]; k <= extent[5]; k++){
                float* ptr = (float*)img->GetScalarPointer(i, j, k);
                if(*ptr == label){
                    count++;
                }
            }
        }
    }
    return count;
}

int main(int argc, char *argv[])
{

    if (argc == 1){
        help(argv[0]);
        return 0;
    }

    string dirname = "";
    string patienfilename = "";
    vector< std::string > patientname;
    string elastixp = "";

    for(int i = 1; i < argc; i++){

        if(string(argv[i]) == "--h" || string(argv[i]) == "--help"){
            help(argv[0]);
            return 0;
        }else if (string(argv[i]) == "-d"){
            if(i + 1 >= argc){
                std::cout << "Patients directory missing"
                          <<std::endl;
                return EXIT_FAILURE;
            }

            dirname = argv[i + 1];


        }else if(string(argv[i]) == "-pfn"){
            patienfilename = argv[i + 1];
        }else if(string(argv[i]) == "-elastixp"){
            elastixp = argv[i + 1];
        }

    }

    patienfilename.insert(0, "/");
    patienfilename.insert(0, dirname);
    ifstream pnfile (patienfilename.c_str());
    cout<<"patient filename: "<<patienfilename<<endl;




    if (pnfile.is_open())
    {
        string line;
        while ( pnfile.good() ){
            getline (pnfile,line);
            patientname.push_back(line);
        }
    }


    vector< vnl_vector< int > > alllabels;
    vector< map< string, int > > alllabelsmap;
    map< string, int > labelsmap;
    vector< vtkSmartPointer< vtkImageData > > alllabelsimages;

    for(unsigned i = 0; i < patientname.size(); i++){

        string annotfilename = dirname;
        annotfilename.append("/");
        annotfilename.append(patientname[i]);
        annotfilename.append("/surf/vtk/");
        annotfilename.append(m_Hemi);
        annotfilename.append(".labelsVertices.txt");


        ifstream annotfile (annotfilename.c_str());

        if (annotfile.is_open())
        {

            vector<int> tempvecannot;

            while ( annotfile.good() ){
                string num;
                getline (annotfile,num);
                if(num.compare("") != 0)
                    tempvecannot.push_back(atoi(num.c_str()));

            }
            annotfile.close();

            vnl_vector<int> vnltemp;
            vnltemp.set_size(tempvecannot.size());
            for(unsigned i=0; i < tempvecannot.size(); i++){
                vnltemp[i] = tempvecannot[i];
            }

            alllabels.push_back(vnltemp);

            tempvecannot.clear();
        }else{
             cout << "Unable to get labeled vertices, use the matlab script to create the labels.";
             return 0;
        }


        string labelannotfilename = dirname;
        labelannotfilename.append("/");
        labelannotfilename.append(patientname[i]);
        labelannotfilename.append("/surf/vtk/");
        labelannotfilename.append(m_Hemi);
        labelannotfilename.append(".labels.txt");


        ifstream labelfile (labelannotfilename.c_str());

        map< string, int > maplabels;

        if (labelfile.is_open()){

            string line;

            while ( labelfile.good() ){

                getline (labelfile,line);
                int pos0 = 0;
                //int pos1 = 0;
                pos0 = line.find(" ");

                if(pos0 != -1){

                    string label = line.substr(0, pos0);

                    int numlabel = atoi(line.substr(pos0+1, line.size() - pos0).c_str());

                    maplabels[label] = numlabel;
                }

            }
            labelfile.close();

        }else{
             cout << "Unable to get annotation file. Use the matlab script to create the labels.";
             return 0;
        }


        alllabelsmap.push_back(maplabels);


        string labelimage = dirname;
        labelimage.append("/");
        labelimage.append(patientname[i]);
        labelimage.append("/surf/vtk/");
        labelimage.append(m_Hemi);
        labelimage.append(".grid.mhd");

        vtkSmartPointer< vtkMetaImageReader > reader = vtkSmartPointer< vtkMetaImageReader >::New();
        reader->SetFileName(labelimage.c_str());
        reader->Update();

        alllabelsimages.push_back(reader->GetOutput());


    }

    for(unsigned i = 0; i < alllabelsmap.size(); i++){
        map< string, int > labelsmapcurrent = alllabelsmap[i];
        map< string, int >::iterator it;
        for(it = labelsmapcurrent.begin(); it != labelsmapcurrent.end(); ++it){
            map< string, int >::iterator fit = labelsmap.find(it->first);
            if(fit != labelsmap.end()){
                if(fit->second != it->second){
                    cout<<"Warning different labels on the structures "<<
                          fit->first<<" "<<fit->second<<" & "<<it->first<<" "<<it->second<<endl;
                }
            }else{
                labelsmap[it->first] = it->second;
            }
        }

    }

    FILE * flabels;
    flabels = fopen ("allLabelsRegistration.txt","w");
    if (flabels!=NULL)
    {

        map< string, int >::iterator it;
        for(it = labelsmap.begin(); it != labelsmap.end(); ++it){

            string currentlabelname = it->first;
            int currentlabel = it->second;

            int minindex = -1;
            int mincount = 99999;
            for(unsigned i = 0; i < alllabelsimages.size(); i++){
                vtkSmartPointer<vtkImageData> img = alllabelsimages[i];
                int count = countLabels(img, currentlabel);

                if(count < mincount){
                    mincount = count;
                    minindex = i;
                }

                string fixedimage = dirname;
                fixedimage.append("/");
                fixedimage.append(patientname[i]);
                fixedimage.append("/surf/vtk/");
                fixedimage.append(currentlabelname);
                fixedimage.append(".maskFixedImage.mhd");
                createMask(img, currentlabel, fixedimage);
            }

            if(minindex != -1){

                fputs (currentlabelname.c_str(),flabels);
                fputs (".",flabels);
                fputs (patientname[minindex].c_str(),flabels);
                fputs (".outputPoints.txt",flabels);

                string elastixexec = elastixp;
                elastixexec.append("/elastix");

                string transformixexec = elastixp;
                transformixexec.append("/transformix");

                string movingimage = dirname;
                movingimage.append("/");
                movingimage.append(patientname[minindex]);
                movingimage.append("/surf/vtk/lh.grid.mhd");

                createGridPoints(alllabelsimages[minindex], currentlabel, alllabelsimages[minindex]->GetNumberOfPoints(), true, "gridPoints.txt");

                string paramfile = elastixp;
                paramfile.append("/");
                paramfile.append("Parameters_BSpline.txt");

                for(unsigned i = 0; i < patientname.size(); i++){

                    if(i != minindex){
                        string fixedimage = dirname;
                        fixedimage.append("/");
                        fixedimage.append(patientname[i]);
                        fixedimage.append("/surf/vtk/lh.grid.mhd");

                        string fixedimageMask = dirname;
                        fixedimageMask.append("/");
                        fixedimageMask.append(patientname[i]);
                        fixedimageMask.append("/surf/vtk/");
                        fixedimageMask.append(currentlabelname);
                        fixedimageMask.append(".maskFixedImage.mhd");

                        string outdir = dirname;
                        outdir.append("/");
                        outdir.append(patientname[i]);
                        outdir.append("/surf/vtk/");


                        string exec = elastixexec;
                        exec.append(" -f ");
                        exec.append(fixedimage);
                        //exec.append(" -fmask ");
                        //exec.append(fixedimageMask);
                        exec.append(" -out ");
                        exec.append(outdir);
                        exec.append(" -m ");
                        exec.append(movingimage);
                        exec.append(" -p ");
                        exec.append(paramfile);

                        cout<<exec<<endl;
                        system(exec.c_str());


                        exec = transformixexec;
                        exec.append(" -def gridPoints.txt -out ");
                        exec.append(outdir);
                        exec.append(" -tp ");
                        exec.append(outdir);
                        exec.append("TransformParameters.0.txt");

                        cout<<exec<<endl;
                        system(exec.c_str());

                        exec = "mv ";
                        exec.append(outdir);
                        exec.append("outputpoints.txt ");
                        exec.append(outdir);
                        exec.append(currentlabelname);
                        exec.append(".");
                        exec.append(patientname[minindex]);
                        exec.append(".outputpoints.txt");
                        cout<<exec<<endl;
                        system(exec.c_str());

                        string filegridpoints = outdir;
                        filegridpoints.append(currentlabelname);
                        filegridpoints.append(".");
                        filegridpoints.append(patientname[minindex]);
                        filegridpoints.append(".gridPoints.txt");
                        createGridPoints(alllabelsimages[minindex], currentlabel, mincount, false, filegridpoints);

                    }
                }

            }

        }

        fclose (flabels);
    }


    return 0;
}


