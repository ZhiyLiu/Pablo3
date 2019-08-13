
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
    cout<<"Test the stereographic projection using a registration."<<endl;
    cout<<"Usage: "<<execname<<" -d <patients dir> -p <patient name> [options]"<<endl;
    cout<<"options:"<<endl;
    cout<<"--h --help show help menu"<<endl;
    cout<<"-p <patient name> add another patient";
    cout<<"-hemi <brain hemisphere lh or rh> default"<<m_Hemi<<endl;
    cout<<"-surf <int> 0 both, 1 top, 2 bottom  surface to be displayed"<<endl;
    cout<<"-min <double> min limit default "<<m_Min<<endl;
    cout<<"-max <double> max limit default "<<m_Max<<endl;
    cout<<"-nump <int> half number of points for one side of the grid the grid default = "<<m_NumP<<", *  2 ="<<m_NumP*2.0<<endl;
    cout<<"-optim <bool> optimize atoms default = "<<m_Optimize<<endl;
    //cout<<"-niter <int> number of iterations for the optimization over the atoms default = "<<m_Niter<<endl;


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
                std::cout << "Patients directory missing"
                          <<std::endl;
                return EXIT_FAILURE;
            }

            dirname = argv[i + 1];


        }else if(string(argv[i]) == "-p"){

            //imagesize[0] = atoi(argv[i+1]);
            patientname.push_back(argv[i + 1]);

        }else if(string(argv[i]) == "-min"){
            m_Min = atof(argv[i + 1]);
        }else if(string(argv[i]) == "-max"){
            m_Max = atof(argv[i + 1]);
        }else if(string(argv[i]) == "-nump"){
            m_NumP = atoi(argv[i + 1]);

        }else if(string(argv[i]) == "-rendersphere"){
            m_RenderSphere = 1;
        }else if(string(argv[i]) == "-startinteractor"){
            m_StartInteractor = 0;
        }else if(string(argv[i]) == "-optim"){
            m_Optimize = atoi(argv[i + 1]);
        }else if(string(argv[i]) == "-hemi"){
            m_Hemi = string(argv[i + 1]);
        }else if(string(argv[i]) == "-interp"){
            m_Inter = atoi(argv[i + 1]);
        }else if(string(argv[i]) == "-surf"){
            m_Surf = atoi(argv[i + 1]);
        }

    }

    if(patientname.size() == 0){
        cout<<"Error Patient name missing!"<<endl;
        return EXIT_FAILURE;
    }

    cout<<"patient: "<<patientname[0]<<endl;
    cout<<"min: "<<m_Min<<endl;
    cout<<"max: "<<m_Max<<endl;
    cout<<"nump: "<<m_NumP<<endl;
    cout<<"hemisphere: "<<m_Hemi<<endl;
    cout<<"Surf: "<<m_Surf<<endl;
    cout<<"interp: "<<m_Inter<<endl;



    vector< vector< float > > allcurvature;
    vnl_vector<float> vectcurvature;
    vector< vnl_vector< int > > alllabels;
    vector< map< string, int > > alllabelsmap;

    vector< vnl_vector<float> > allthickness;

    for(unsigned i = 0; i < patientname.size(); i++){
        string curvfilename = dirname;
        curvfilename.append("/");
        curvfilename.append(patientname[i]);
        curvfilename.append("/surf/vtk/");
        curvfilename.append(m_Hemi);
        curvfilename.append(".curv.asc");

        string line;
        ifstream myfile (curvfilename.c_str());

        if (myfile.is_open())
        {

            vector<float> tempveccurvature;

            while ( myfile.good() ){
                getline (myfile,line);
                int pos0 = 0;
                pos0 = line.find_last_of(" ");

                if(pos0 != -1){
                    string num = line.substr(pos0, line.size() - pos0);


                    tempveccurvature.push_back(atof(num.c_str()));
                }



            }
            myfile.close();

            if(vectcurvature.size() == 0){
                vectcurvature.set_size(tempveccurvature.size());
                vectcurvature.fill(0);
            }

            for(unsigned i = 0; i < vectcurvature.size(); i++){
                vectcurvature[i] += tempveccurvature[i];
            }

            allcurvature.push_back(tempveccurvature);

            tempveccurvature.clear();
        }else{
             cout << "Unable to load the curvatures from "<<curvfilename<<endl;
             return 0;
        }



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



        string thickfilename = dirname;
        thickfilename.append("/");
        thickfilename.append(patientname[i]);
        thickfilename.append("/surf/vtk/");
        thickfilename.append(m_Hemi);
        thickfilename.append(".thickness.asc");

        ifstream thickfile (thickfilename.c_str());

        if (thickfile.is_open())
        {

            vnl_vector<float> vectthickness;
            vector<float> tempvecthick;

            while ( thickfile.good() ){
                getline (thickfile,line);
                int pos0 = 0;
                pos0 = line.find_last_of(" ");

                if(pos0 != -1){
                    string num = line.substr(pos0, line.size() - pos0);


                    tempvecthick.push_back(atof(num.c_str()));
                }



            }
            thickfile.close();

            if(vectthickness.size() == 0){
                vectthickness.set_size(tempvecthick.size());
                vectthickness.fill(0);
            }

            for(unsigned i = 0; i < tempvecthick.size(); i++){
                vectthickness[i] = tempvecthick[i];
            }

            allthickness.push_back(vectthickness);

            vectthickness.clear();
            tempvecthick.clear();
        }else{
             cout << "Unable to label thickness file from freesurfer.";
             return 0;
        }
    }


    vectcurvature = vectcurvature / (patientname.size());


    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1,1,1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    vtkSmartPointer<vtkTransform> transformaxes = vtkSmartPointer<vtkTransform>::New();
    transformaxes->Translate(0,0, 100.0);
    transformaxes->Scale(10,10,10);
    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
      // The axes are positioned with a user transform
      axes->SetUserTransform(transformaxes);
      axes->SetAxisLabels(0);

      // properties of the axes labels can be set as follows
      // this sets the x axis label to red
      // axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1,0,0);

      // the actual text of the axis label can be changed:
      // axes->SetXAxisLabelText("test");

      //renderer->AddActor(axes);





    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();

    string origfilename = dirname;
    origfilename.append("/");
    origfilename.append(patientname[0]);
    origfilename.append("/surf/vtk/");
    origfilename.append(m_Hemi);
    origfilename.append(".orig.vtk");

    reader->SetFileName(origfilename.c_str());
    reader->Update();
    vtkPolyData* orig = reader->GetOutput();



    string pialfilename = dirname;
    pialfilename.append("/");
    pialfilename.append(patientname[0]);
    pialfilename.append("/surf/vtk/");
    pialfilename.append(m_Hemi);
    pialfilename.append(".pial.vtk");

    vtkSmartPointer<vtkPolyDataReader> reader2 =  vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(pialfilename.c_str());
    reader2->Update();
    vtkPolyData* pial = reader2->GetOutput();

    vtkSmartPointer<vtkPolyDataNormals> polynormals = vtkSmartPointer<vtkPolyDataNormals>::New();
    polynormals->SetInput(pial);
    polynormals->ComputePointNormalsOn();
    polynormals->Update();
    vtkPolyData* brainnormals = polynormals->GetOutput();

    vtkSmartPointer<vtkPolyDataReader> reader3 =  vtkSmartPointer<vtkPolyDataReader>::New();

    string spherefilename = dirname;
    spherefilename.append("/");
    spherefilename.append(patientname[0]);
    spherefilename.append("/surf/vtk/");
    spherefilename.append(m_Hemi);
    spherefilename.append(".sphere.vtk");

    reader3->SetFileName(spherefilename.c_str());
    reader3->Update();
    vtkPolyData* sphere = reader3->GetOutput();


    //vtkDataArray* labelarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
    //double* curvrange = curv->GetOutput()->GetPointData()->GetScalars()->GetRange();

    vnl_vector<double> northpole;
    double angleplane = 0;
    vector< vnl_vector< double > > sampleborderpoints;

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();

    vnl_vector<double> north(3);
    north.fill(0);
    north[2] = 1;

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFiltercircle =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    vtkSmartPointer<vtkPolyData> pointsprojectedcircle = 0;
    vtkSmartPointer<vtkPolyData> bestcircle = 0;
    vtkSmartPointer<vtkPolyData> bestcircletangent = 0;



    vtkSmartPointer<vtkColorTransferFunction> colorannot = vtkSmartPointer<vtkColorTransferFunction>::New();
    colorannot->SetColorSpaceToHSV();
    vtkDataArray* labelarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);

    for(unsigned i = 0; i < 1; i++){
        vnl_vector<int> vectlabel = alllabels[i];
        map< string, int > maplabels = alllabelsmap[i];


        double step = 1.0/((double)maplabels.size());
        double hue = 0;

        map< string, int >::iterator it;

        for(unsigned j = 0; j < vectlabel.size(); j++){
            double lab[1];
            lab[0] = vectlabel[j];
            labelarray->InsertTuple(j, lab);
        }

        for ( it=maplabels.begin() ; it != maplabels.end(); it++ ){

            int label = (*it).second;
            colorannot->AddHSVPoint(label, hue, 0.6, 1);
            hue += step;

        }

        int corpuslabel = maplabels["corpuscallosum"];

        vtkSmartPointer< vtkBestCircle > findbestcircle = vtkSmartPointer< vtkBestCircle >::New();
        findbestcircle->SetInput(sphere);
        findbestcircle->SetLabels(vectlabel);
        findbestcircle->SetLabel(corpuslabel);
        findbestcircle->Update();



        northpole = findbestcircle->GetCircleCenter();
        angleplane = -findbestcircle->GetAnglePlane()*180.0/PI;
        sampleborderpoints = findbestcircle->GetSampleBorderPoints();
        vector< vnl_vector< double > > sampleborderpoints = findbestcircle->GetSampleBorderPoints();
        vtkSmartPointer< vtkInterpolateCurve > interpolatecurve = vtkSmartPointer< vtkInterpolateCurve >::New();
        interpolatecurve->SetCurvePoints(sampleborderpoints);
        interpolatecurve->SetCyclic(true);
        interpolatecurve->SetInterpolationLevel(3);
        interpolatecurve->Update();
        vtkSmartPointer<vtkPolyDataMapper> curvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        curvemapper->SetInput(interpolatecurve->GetOutput());
        vtkSmartPointer<vtkActor>  vtkactorcurve = vtkActor::New();
        vtkactorcurve->SetMapper(curvemapper);

        vtkactorcurve->GetProperty()->SetLineWidth(4);
        vtkactorcurve->GetProperty()->SetColor(0,1,1);
        //renderer->AddActor(vtkactorcurve);



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



    }

    colorannot->Build();


    vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
    color->SetColorSpaceToRGB();

    //double* curvrange = curv->GetOutput()->GetPointData()->GetScalars()->GetRange();

    color->AddRGBPoint(-0.2,1,0,0);
    color->AddRGBPoint(0.2,0,1,0);

    color->Build();

    vtkDataArray* curvaturearray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);

    for(unsigned i = 0; i < vectcurvature.size();i++){
        double val[1];
        val[0] = vectcurvature[i];
        curvaturearray->InsertTuple(i, val);
    }


    vtkSmartPointer<vtkTransform> transform1 = vtkSmartPointer<vtkTransform>::New();
    transform1->RotateWXYZ(angleplane, north[0], north[1], north[2]);
    //cout<<"angle "<<angleplane<<endl;
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter1 =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter1->SetTransform(transform1);
    transformFilter1->SetInputConnection(sphere->GetProducerPort());
    transformFilter1->Update();
    sphere = transformFilter1->GetOutput();

    /*vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter2 =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter2->SetTransform(transform1);
    transformFilter2->SetInputConnection(pointsprojectedcircle->GetProducerPort());
    transformFilter2->Update();
    pointsprojectedcircle = transformFilter2->GetOutput();


    vtkSmartPointer<vtkPolyDataMapper> pointsprojectedcirclemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointsprojectedcirclemapper->SetInput(pointsprojectedcircle);
    vtkSmartPointer<vtkActor>  vtkactorpointsprojectedcircle = vtkActor::New();
    vtkactorpointsprojectedcircle->SetMapper(pointsprojectedcirclemapper);
    vtkactorpointsprojectedcircle->GetProperty()->SetColor(0,1,1);
    vtkactorpointsprojectedcircle->GetProperty()->SetPointSize(6);
    renderer->AddActor(vtkactorpointsprojectedcircle);


    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilterbestcircle =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilterbestcircle->SetTransform(transform);
    transformFilterbestcircle->SetInput(bestcircle);
    bestcircle = transformFilterbestcircle->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> bestcirclemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    bestcirclemapper->SetInput(bestcircle);
    vtkSmartPointer<vtkActor>  vtkactorbestcircle = vtkActor::New();
    vtkactorbestcircle->SetMapper(bestcirclemapper);
    vtkactorbestcircle->GetProperty()->SetColor(1,1,0);
    vtkactorbestcircle->GetProperty()->SetPointSize(6);
    renderer->AddActor(vtkactorbestcircle);

    //renderer->AddActor(vtkactorcircle);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilterbestcircletangent =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilterbestcircletangent->SetTransform(transform);
    transformFilterbestcircletangent->SetInput(bestcircletangent);
    bestcircletangent = transformFilterbestcircletangent->GetOutput();

    {
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilterbestcircletangent =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilterbestcircletangent->SetTransform(transform1);
        transformFilterbestcircletangent->SetInput(bestcircletangent);
        bestcircletangent = transformFilterbestcircletangent->GetOutput();
    }


    vtkSmartPointer<vtkPolyDataMapper> bestcircletangentmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    bestcircletangentmapper->SetInput(bestcircletangent);
    vtkSmartPointer<vtkActor>  vtkactorbestcircletangent = vtkActor::New();
    vtkactorbestcircletangent->SetMapper(bestcircletangentmapper);
    vtkactorbestcircletangent->GetProperty()->SetColor(1,0,1);
    vtkactorbestcircletangent->GetProperty()->SetLineWidth(5);
    vtkactorbestcircletangent->GetProperty()->SetPointSize(6);
    renderer->AddActor(vtkactorbestcircletangent);*/


    //sphere->GetPointData()->SetScalars(curvaturearray);
    sphere->GetPointData()->SetScalars(labelarray);
    vtkSmartPointer<vtkPolyDataMapper> spheremapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    spheremapper->SetInput(sphere);
    vtkSmartPointer<vtkActor>  vtkactorsphere = vtkActor::New();
    vtkactorsphere->SetMapper(spheremapper);
    spheremapper->SetLookupTable(colorannot);


    vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
    locator->SetDataSet(sphere);
    locator->AutomaticOn();
    locator->SetNumberOfPointsPerBucket(3);
    locator->BuildLocator();


    vtkDataArray* normalarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
    normalarray->SetNumberOfComponents(3);


    string transformpoints = dirname;
    transformpoints.append("/");
    transformpoints.append(patientname[0]);
    transformpoints.append("/surf/vtk/");
    transformpoints.append("outputpoints.txt");

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
    /*vtkSmartPointer<vtkPolyDataMapper> srepmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    srepmapper->SetInput(srepfig);
    vtkSmartPointer<vtkActor>  srepactor = vtkActor::New();
    srepactor->SetMapper(srepmapper);
    srepmapper->SetLookupTable(colorannot);
    renderer->AddActor(srepactor);*/



    if(m_Optimize){

        /*vtkSmartPointer< vtkOptimizeCortexAtomsCurv > optimizeatomscurv = vtkSmartPointer< vtkOptimizeCortexAtomsCurv >::New();
        optimizeatomscurv->SetInput(srepfig);
        optimizeatomscurv->SetInputSphere(sphere);
        optimizeatomscurv->SetInputCurvatureArray(curvaturearray);
        if(!m_RenderSphere){
            //optimizeatomscurv->SetInputCortexSurf(pial, orig);
        }
        optimizeatomscurv->Update();
        srepfig = optimizeatomscurv->GetOutput();*/

        /*vtkSmartPointer< vtkOptimizeCortexEndAtomsCorpus > optimizeendatoms = vtkSmartPointer< vtkOptimizeCortexEndAtomsCorpus >::New();
        optimizeendatoms->SetInput(srepfig);
        optimizeendatoms->SetInputSphere(sphere);
        optimizeendatoms->SetInputSamplePoints(sampleborderpoints);

        int limit = 2;

        if(!m_RenderSphere && limit == 0){
            optimizeendatoms->SetInputCortexSurf(pial, orig);
        }
        optimizeendatoms->Update();
        srepfig = optimizeendatoms->GetOutput();

        //vtkSmartPointer<vtkPolyDataMapper> corpuscontourmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        //corpuscontourmapper->SetInput(optimizeendatoms->GetOuputCorpusContour());
        //vtkSmartPointer<vtkActor>  vtkactorcorpus = vtkActor::New();
        //vtkactorcorpus->SetMapper(corpuscontourmapper);
        //vtkactorcorpus->GetProperty()->SetLineWidth(3);
        //renderer->AddActor(vtkactorcorpus);*/


        int limit = 1;
        for(unsigned i = 0; i < limit; i++){

            vtkSmartPointer< vtkOptimizeCortexAtomsDis > optimizeatoms = vtkSmartPointer< vtkOptimizeCortexAtomsDis >::New();
            optimizeatoms->SetInput(srepfig);
            optimizeatoms->SetInputSphere(sphere);
            optimizeatoms->SetLabelArray(labelarray);

            if(!m_RenderSphere && i + 1 == limit){
                optimizeatoms->SetInputCortexSurf(pial, orig);
            }
            optimizeatoms->Update();
            srepfig = optimizeatoms->GetOutput();
        }

    }

    cout<<"npoints"<<(int)srepfig->GetNumberOfPoints()<<endl;
    cout<<"ncells"<<(int)srepfig->GetNumberOfCells()<<endl;
    //srepfig->SetGridTopolgyIds(pointsIds);



    vtkSmartPointer<vtkPolyDataMapper> origmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    origmapper->SetInput(orig);
    orig->GetPointData()->SetScalars(labelarray);
    origmapper->SetLookupTable(colorannot);

    vtkSmartPointer<vtkActor>  vtkactororig = vtkActor::New();
    vtkactororig->SetMapper(origmapper);



    //pial->GetPointData()->SetScalars(curvaturearray);
    vtkSmartPointer<vtkPolyDataMapper> pialmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pialmapper->SetInput(pial);
    pial->GetPointData()->SetScalars(labelarray);
    pialmapper->SetLookupTable(colorannot);
    vtkSmartPointer<vtkActor>  vtkactorpial = vtkActor::New();
    vtkactorpial->SetMapper(pialmapper);
    //vtkactorpial->GetProperty()->SetOpacity(0.6);
    //vtkactorpial->GetProperty()->SetColor(0.2,0.5,0.8);
    //pialmapper->SetLookupTable(color);

    double bounds2[6];
    pial->GetBounds(bounds2);
    vnl_vector<double> vnlbounds(bounds2, 6);
    cout<<"bounds "<<vnlbounds<<endl;


    vtkSmartPointer< vtkSRepVisuPrimitives > visuprimitives = vtkSmartPointer< vtkSRepVisuPrimitives >::New();
    visuprimitives->SetInput(srepfig);
    visuprimitives->Update();

    vtkActorCollection* actorcollection = visuprimitives->GetOuputActors();

    for(unsigned i = 0; i < 1; i++){

        vtkActor* actor = (vtkActor*) actorcollection->GetItemAsObject(i);

        renderer->AddActor(actor);
    }

    int interpolationlevel = 0;
    if(!m_RenderSphere){
        interpolationlevel = m_Inter;
    }


    vtkSmartPointer< vtkSRepInterpolateMedialSheet > medialsheetinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSheet >::New();
    medialsheetinterpolator->SetInput(srepfig);
    medialsheetinterpolator->SetInterpolationLevel(interpolationlevel);
    medialsheetinterpolator->Update();


    vtkSmartPointer< vtkPolyData > medialsheet = medialsheetinterpolator->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> medialsheetmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    medialsheetmapper->SetInputConnection(medialsheet->GetProducerPort());
    vtkSmartPointer<vtkActor>  medialsheetactor = vtkActor::New();
    medialsheetactor->SetMapper(medialsheetmapper);
    medialsheetactor->GetProperty()->SetPointSize(10);
    medialsheetactor->GetProperty()->SetLineWidth(3);
    medialsheetactor->GetProperty()->SetColor(0.2,0.6,1);
    medialsheetactor->GetProperty()->SetRepresentationToWireframe();


    vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();
    medialspokesinterpolator->SetInput(srepfig);
    medialspokesinterpolator->SetInterpolationLevel(interpolationlevel);
    if(m_Surf == 1){
        medialspokesinterpolator->SetSpokeType(vtkSRep::TOP_SPOKE);
    }else if(m_Surf == 2){
        medialspokesinterpolator->SetSpokeType(vtkSRep::BOTTOM_SPOKE);
    }
    medialspokesinterpolator->Update();

    vtkPolyData* interpolatedmedialspokes = medialspokesinterpolator->GetOutput();
    /*vtkSmartPointer<vtkPolyDataMapper> srepmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    srepmapper->SetInput(medialspokesinterpolator->GetSRepOutput());
    vtkSmartPointer<vtkActor>  srepactor = vtkActor::New();
    srepactor->SetMapper(srepmapper);
    srepmapper->SetLookupTable(colorannot);
    renderer->AddActor(srepactor);*/


    vtkSmartPointer<vtkPolyDataMapper> medialspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    medialspokesmapper->SetInputConnection(interpolatedmedialspokes->GetProducerPort());
    medialspokesmapper->SetLookupTable(colorannot);
    vtkSmartPointer<vtkActor>  medialspokesactor = vtkActor::New();
    medialspokesactor->SetMapper(medialspokesmapper);
    medialspokesactor->GetProperty()->SetColor(0,1.0,0);
    medialspokesactor->GetProperty()->SetLineWidth(5);


    //medialspokesactor->GetProperty()->SetOpacity(0.8);

    renderer->AddActor(medialsheetactor);

    if(m_RenderSphere){
        //vtkactorsphere->GetProperty()->SetOpacity(0.5);
        renderer->AddActor(vtkactorsphere);
    }else{
        //renderer->AddActor(medialspokesactor);
        if(m_Surf == 1){
            //renderer->AddActor(vtkactorpial);
            //vtkactorpial->GetProperty()->SetOpacity(0.5);
        }else if(m_Surf == 2){
            //renderer->AddActor(vtkactororig);
            //vtkactororig->GetProperty()->SetOpacity(0.5);
        }
    }


    /*vtkSmartPointer< vtkSRepToGridLayout > sreptogrid = vtkSmartPointer< vtkSRepToGridLayout >::New();
    sreptogrid->SetInput(srepfig);
    sreptogrid->Update();

    vtkSmartPointer<vtkPolyDataMapper> srepmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSRep* srepgrid = sreptogrid->GetOutput();
    srepmapper->SetInput(srepgrid);
    srepgrid->GetPointData()->SetScalars(labelarray);
    vtkSmartPointer<vtkActor>  srepactor = vtkActor::New();
    srepactor->SetMapper(srepmapper);
    srepmapper->SetLookupTable(colorannot);
    srepactor->GetProperty()->SetRepresentationToWireframe();
    srepactor->GetProperty()->SetLineWidth(4);
    renderer->AddActor(srepactor);*/



    if(m_StartInteractor)
        renderWindowInteractor->Start();




    return 0;
}


