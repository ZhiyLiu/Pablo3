
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


double m_Min = 0.04;
double m_Max = 1;
int m_NumP = 40;
int m_Inter = 3;
bool m_RenderSphere = 0;
bool m_StartInteractor = 1;
bool m_Optimize = 1;
string m_Hemi = "lh";
int m_Surf = 0;
string m_FileAllLabels = "allLabelsRegistration.txt";

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
    cout<<"Creates an s-rep per structure for all the subjects"<<endl;
    cout<<"Usage: "<<execname<<" -d <patients dir> -pfn <patient file name> [options]"<<endl;
    cout<<"options:"<<endl;
    cout<<"--h --help show help menu"<<endl;
    cout<<"-pfn <patient filename> list of file with name of patients <one per line>"<<endl;
    cout<<"-hemi <brain hemisphere lh or rh> default"<<m_Hemi<<endl;
    cout<<"-min <double> min limit default "<<m_Min<<endl;
    cout<<"-max <double> max limit default "<<m_Max<<endl;
    cout<<"-nump <int> half number of points for one side of the grid the grid default = "<<m_NumP<<", *  2 ="<<m_NumP*2.0<<endl;
    cout<<"-labels <file with the labels for the structures default: "<<m_FileAllLabels<<"(see testRegistration)"<<endl;


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
        }else if(string(argv[i]) == "-min"){
            m_Min = atof(argv[i + 1]);
        }else if(string(argv[i]) == "-max"){
            m_Max = atof(argv[i + 1]);
        }else if(string(argv[i]) == "-nump"){
            m_NumP = atoi(argv[i + 1]);

        }else if(string(argv[i]) == "-hemi"){
            m_Hemi = string(argv[i + 1]);
        }else if(string(argv[i]) == "-labels"){
            m_FileAllLabels = string(argv[i + 1]);
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

    for(unsigned i = 0; i < patientname.size(); i++){

        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();

        string origfilename = dirname;
        origfilename.append("/");
        origfilename.append(patientname[i]);
        origfilename.append("/surf/vtk/");
        origfilename.append(m_Hemi);
        origfilename.append(".orig.vtk");

        reader->SetFileName(origfilename.c_str());
        reader->Update();
        vtkPolyData* orig = reader->GetOutput();


        string pialfilename = dirname;
        pialfilename.append("/");
        pialfilename.append(patientname[i]);
        pialfilename.append("/surf/vtk/");
        pialfilename.append(m_Hemi);
        pialfilename.append(".pial.vtk");

        vtkSmartPointer<vtkPolyDataReader> reader2 =  vtkSmartPointer<vtkPolyDataReader>::New();
        reader2->SetFileName(pialfilename.c_str());
        reader2->Update();
        vtkPolyData* pial = reader2->GetOutput();

        vtkSmartPointer<vtkPolyDataReader> reader3 =  vtkSmartPointer<vtkPolyDataReader>::New();
        string spherefilename = dirname;
        spherefilename.append("/");
        spherefilename.append(patientname[i]);
        spherefilename.append("/surf/vtk/");
        spherefilename.append(m_Hemi);
        spherefilename.append(".sphere.vtk");

        reader3->SetFileName(spherefilename.c_str());
        reader3->Update();
        vtkPolyData* sphere = reader3->GetOutput();

        vnl_vector<double> northpole;
        double angleplane = 0;

        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();

        vnl_vector<double> north(3);
        north.fill(0);
        north[2] = 1;

        vtkSmartPointer<vtkTransformPolyDataFilter> transformFiltercircle =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        vtkSmartPointer<vtkPolyData> pointsprojectedcircle = 0;
        vtkSmartPointer<vtkPolyData> bestcircle = 0;
        vtkSmartPointer<vtkPolyData> bestcircletangent = 0;


        vnl_vector<int> vectlabel = alllabels[i];
        int corpuslabel = labelsmap["corpuscallosum"];

        vtkSmartPointer< vtkBestCircle > findbestcircle = vtkSmartPointer< vtkBestCircle >::New();
        findbestcircle->SetInput(sphere);
        findbestcircle->SetLabels(vectlabel);
        findbestcircle->SetLabel(corpuslabel);
        findbestcircle->Update();


        northpole = findbestcircle->GetCircleCenter();
        angleplane = -findbestcircle->GetAnglePlane()*180.0/PI;


        vnl_vector<double> cross = vnl_cross_3d(northpole, north);
        double angle = acos(dot_product(northpole, north))*180.0/PI;

        transform->RotateWXYZ(angle, cross[0], cross[1], cross[2]);


        transformFilter->SetTransform(transform);
        transformFilter->SetInputConnection(sphere->GetProducerPort());
        transformFilter->Update();

        sphere = transformFilter->GetOutput();

        transformFiltercircle->SetTransform(transform);
        transformFiltercircle->SetInputConnection(findbestcircle->GetOutput()->GetProducerPort());
        transformFiltercircle->Update();

        pointsprojectedcircle = transformFiltercircle->GetOutput();
        //pointsprojectedcircle = findbestcircle->GetOutput();

        bestcircle = findbestcircle->GetOutputCircle();
        bestcircletangent = findbestcircle->GetOutputTangent();


        vtkSmartPointer<vtkTransform> transform1 = vtkSmartPointer<vtkTransform>::New();
        transform1->RotateWXYZ(angleplane, north[0], north[1], north[2]);
        //cout<<"angle "<<angleplane<<endl;
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter1 =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter1->SetTransform(transform1);
        transformFilter1->SetInputConnection(sphere->GetProducerPort());
        transformFilter1->Update();
        sphere = transformFilter1->GetOutput();



        ifstream alllabelsfile (m_FileAllLabels.c_str());

        if (alllabelsfile.is_open()){

            string labelfile;

            while ( alllabelsfile.good() ){

                getline (alllabelsfile, labelfile);


                string transformpoints = dirname;
                transformpoints.append("/");
                transformpoints.append(patientname[i]);
                transformpoints.append("/surf/vtk/");
                transformpoints.append(labelfile);

                string line;
                ifstream myfile (transformpoints.c_str());
                vector< vnl_vector<double> > registeredpoints;

                if (myfile.is_open()){


                    while ( myfile.good() ){
                        getline (myfile,line);
                        int pos0 = 0, pos1 = 0;
                        string tof = "OutputPoint = [ ";
                        pos0 = line.find(tof) + tof.size();
                        pos1 = line.find(" ]", pos0);

                        if(pos0 != -1 && pos1 != -1){
                            string num = line.substr(pos0, pos1 - pos0);
                            string temp = num.substr(0, num.find(" "));
                            double x = atof(temp.c_str());
                            temp = num.substr(num.find(" "), num.size());
                            double y = atof(temp.c_str());

                            vnl_vector<double> p(3);
                            p.fill(0);
                            p[0] = x;
                            p[1] = y;
                            registeredpoints.push_back(p);
                        }

                    }
                    myfile.close();
                }


                vtkSmartPointer< vtkOptimizeSRepSampling > srepsampling = vtkSmartPointer< vtkOptimizeSRepSampling >::New();
                srepsampling->SetMax(m_Max);
                srepsampling->SetMin(m_Min);
                srepsampling->SetNumP(m_NumP);
                srepsampling->SetRegisteredPoints(registeredpoints);
                srepsampling->Update();
                vtkSmartPointer<vtkSRep> srepfig = srepsampling->GetOutput();

                vtkSmartPointer< vtkOptimizeCortexAtomsDis > optimizeatoms = vtkSmartPointer< vtkOptimizeCortexAtomsDis >::New();
                optimizeatoms->SetInput(srepfig);
                optimizeatoms->SetInputSphere(sphere);
                optimizeatoms->SetInputCortexSurf(pial, orig);
                optimizeatoms->Update();
                srepfig = optimizeatoms->GetOutput();


                cout<<"npoints"<<(int)srepfig->GetNumberOfPoints()<<endl;
                cout<<"ncells"<<(int)srepfig->GetNumberOfCells()<<endl;


                M3DQuadFigure* cortex = NewQuadFigure(srepfig->GetNumColumns(), srepfig->GetNumRows());

                for(unsigned l = 0; l < srepfig->GetNumColumns(); l++){
                     for(unsigned m = 0; m < srepfig->GetNumRows(); m++){

                         double point[3];
                         vtkIdType srepspokeid = l*srepfig->GetNumColumns() + m;
                         srepfig->GetPoint(srepspokeid, point);

                         M3DQuadPrimitive* prim = dynamic_cast<M3DQuadPrimitive*>(cortex->getPrimitivePtr(l, m));
                         prim->setX(point[0], point[1], point[2]);

                         vtkSRep::VNLType vnlvect0 = srepfig->GetSpoke(srepspokeid, vtkSRep::TOP_SPOKE);
                         Vector3D vect0 = prim->getU0();
                         vect0.setX(vnlvect0[0]);
                         vect0.setY(vnlvect0[1]);
                         vect0.setZ(vnlvect0[2]);
                         prim->setU0(vect0);
                         prim->setR(srepfig->GetSpokeRadius(srepspokeid, vtkSRep::TOP_SPOKE));


                         vtkSRep::VNLType vnlvect1 = srepfig->GetSpoke(srepspokeid, vtkSRep::BOTTOM_SPOKE);
                         Vector3D vect1 = prim->getU1();
                         vect1.setX(vnlvect1[0]);
                         vect1.setY(vnlvect1[1]);
                         vect1.setZ(vnlvect1[2]);
                         prim->setU1(vect1);
                         prim->setR(srepfig->GetSpokeRadius(srepspokeid, vtkSRep::BOTTOM_SPOKE));

                         Vector3D vectend = prim->getUEnd();
                         vectend.setX(0);
                         vectend.setY(0);
                         vectend.setZ(0);
                         prim->setUEnd(vectend);

                     }
                 }

                string m3dfilename = dirname;
                m3dfilename.append("/");
                m3dfilename.append(patientname[i]);
                m3dfilename.append("/surf/vtk/");
                m3dfilename.append(m_Hemi);
                m3dfilename.append(".");
                m3dfilename.append(labelfile.substr(0, labelfile.find_first_of(".")));
                m3dfilename.append(".medial.m3d");

                SaveQuadFigure(cortex, m3dfilename.c_str());

            }
            alllabelsfile.close();
        }
    }

    return 0;
}


