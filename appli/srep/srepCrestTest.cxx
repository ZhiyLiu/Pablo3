#include <vtkSmartPointer.h>

#include <vtksrep.h>
#include <vtksrepinterpolatemedialsheet.h>
#include <vtksrepvisuprimitives.h>

#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkLine.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#include <vtksrepinterpolatemedialcrestcurve.h>
#include <vtksrepinterpolatemedialspokes.h>
#include <vtksrepinterpolatemedialspokeshermite.h>
#include <vtksrepinterpolatecrestspokes.h>
#include <vtksrepinterpolatecrestspokesquartic.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkCleanPolyData.h>

#include "P3DControl.h"
#include "ControlParms.h"

#include <vtkInterpolateDataSetAttributes.h>

ControlParms * globalControl;	// Read the user's preferences file
int globalVerbosity;			// Current verbosity level of Pablo

using namespace std;

void help(char* execname){
    cout<<"Srep test using vtk"<<endl;
    //cout<<"Usage: "<<execname<<" -d <patients dir> -p <patient name> [options]"<<endl;
    //cout<<"options:"<<endl;
    cout<<"--h --help show help menu"<<endl;
    cout<<"-m <model name> ";

}

M3DQuadFigure* GetQuadFigure(const char* figfilename){
    globalControl = new ControlParms(NULL, -1, false);	// Ignore user's preferences
    globalVerbosity = 0;
    globalControl->setDefault(OutputVerbosity, globalVerbosity);
    globalControl->setDefault(ReorderModels, 0);
    globalControl->setDefault(SmoothImages, false);
    globalControl->setDefault(ConvertImages, false);
    globalControl->setDefault(ByteOrder, 1);
    globalControl->setDefault(CompressImages, true);
    globalControl->setDefault(ShowLandmarks, true);

    P3DControl* p3d = new P3DControl(10);
    p3d->read(figfilename, false);
    M3DObject* m3dobject = p3d->getObjectPtr();
    M3DFigure * figure = m3dobject->getFigurePtr( 0 ) ;
   return dynamic_cast<M3DQuadFigure*>( figure );
}

int main(int argc, char *argv[])
{


    if (argc < 1){
        help(argv[0]);
        return 0;
    }


    vector< std::string > modelname;

    for(int i = 1; i < argc; i++){

        if(string(argv[i]) == "--h" || string(argv[i]) == "--help"){
            help(argv[0]);
            return 0;
        }else if(string(argv[i]) == "-m"){
            modelname.push_back(argv[i + 1]);
        }

    }

         vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
         renderer->SetBackground(1,1,1);
         vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
         renderWindow->AddRenderer(renderer);
         vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
         renderWindowInteractor->SetRenderWindow(renderWindow);


         for(unsigned modi = 0; modi < modelname.size(); modi++){

             M3DQuadFigure* quadfig = GetQuadFigure(modelname[modi].c_str());



             vtkSmartPointer<vtkSRep> srepfig = vtkSmartPointer<vtkSRep>::New();
             vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
             vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

             vtkSRep::VectorSRepIdsType pointsIds;

             vtkSRep::RadiusVectorType allradius;
             vtkSRep::SpokesVectorType allspokes;

             for(int u = 0; u < quadfig->getRowCount(); u++){
                 pointsIds.push_back(vtkSRep::VectorIdsType());
                 for(int v = 0; v < quadfig->getColumnCount(); v++){

                     M3DQuadPrimitive* prim0 = dynamic_cast<M3DQuadPrimitive*>(quadfig->getPrimitivePtr(u, v));

                     Vector3D x = prim0->getX();
                     Vector3D u0 = prim0->getU0();
                     Vector3D u1 = prim0->getU1();

                     vtkSRep::VectorVNLType vnlspokes;
                     vtkSRep::VNLType s(3);
                     s[0] = u0.getX();
                     s[1] = u0.getY();
                     s[2] = u0.getZ();
                     vnlspokes.push_back(s);

                     s[0] = u1.getX();
                     s[1] = u1.getY();
                     s[2] = u1.getZ();
                     vnlspokes.push_back(s);

                     vtkSRep::VectorDoubleType radius;
                     radius.push_back(prim0->getR0());
                     radius.push_back(prim0->getR1());

                     if(u == 0 || u == quadfig->getRowCount() - 1 || v == 0 || v == quadfig->getColumnCount() - 1){

                         M3DQuadEndPrimitive* prim0 = dynamic_cast<M3DQuadEndPrimitive*>(quadfig->getPrimitivePtr(u, v));
                         Vector3D uend = prim0->getUEnd();

                         s[0] = uend.getX();
                         s[1] = uend.getY();
                         s[2] = uend.getZ();

                         vnlspokes.push_back(s);

                         radius.push_back(prim0->getREnd());

                     }


                     pointsIds[u].push_back(hubpos->InsertNextPoint(x.getX(), x.getY(), x.getZ()));



                     allspokes.push_back(vnlspokes);
                     allradius.push_back(radius);
                 }
             }


             for(unsigned i = 0; i < pointsIds.size() - 1; i++){
                  for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

                      vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                      quad->GetPointIds()->SetId(0, pointsIds[i][j]);
                      quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
                      quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
                      quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

                      //quad->Print(cout);

                      cellarray->InsertNextCell(quad);

                  }
              }

             srepfig->SetPoints(hubpos);
             srepfig->SetPolys(cellarray);
             srepfig->SetAllSpokes(allspokes);
             srepfig->SetAllSpokesRadius(allradius);
             //srepfig->SetGridTopolgyIds(pointsIds);


             int interpolationlevel = 3;


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
             medialsheetactor->GetProperty()->SetColor(0.2,0.6,1);
             renderer->AddActor(medialsheetactor);


            //vtkSmartPointer< vtkSRepInterpolateMedialSpokes > medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokes >::New();
            vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();

            medialspokesinterpolator->SetInput(srepfig);
            medialspokesinterpolator->SetInterpolationLevel(interpolationlevel);
            //medialspokesinterpolator->SetAtomId(4);
            //medialspokesinterpolator->SetGamma_u(0.5);
            //medialspokesinterpolator->SetGamma_v(0.5);
            medialspokesinterpolator->Update();
            vtkPolyData* interpolatedmedialspokes = medialspokesinterpolator->GetOutput();

            vtkSmartPointer<vtkPolyDataMapper> medialspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            medialspokesmapper->SetInputConnection(interpolatedmedialspokes->GetProducerPort());
            vtkSmartPointer<vtkActor>  medialspokesactor = vtkActor::New();
            medialspokesactor->SetMapper(medialspokesmapper);
            medialspokesactor->GetProperty()->SetColor(1,0.9,0.1);
            medialspokesactor->GetProperty()->SetLineWidth(5);
            //renderer->AddActor(medialspokesactor);


             vtkSmartPointer<vtkSRepInterpolateMedialCrestCurve> curveinterpolation = vtkSmartPointer<vtkSRepInterpolateMedialCrestCurve>::New();
             curveinterpolation->SetInput(srepfig);
             curveinterpolation->SetInterpolationLevel(interpolationlevel);
             curveinterpolation->Update();
             vtkPolyData* medialcrestcurve = curveinterpolation->GetOutput();

             vtkSmartPointer<vtkPolyDataMapper> medialcrestcurvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
             medialcrestcurvemapper->SetInputConnection(medialcrestcurve->GetProducerPort());
             vtkSmartPointer<vtkActor>  medialsheetcrestcurveactor = vtkActor::New();
             medialsheetcrestcurveactor->SetMapper(medialcrestcurvemapper);
             medialsheetcrestcurveactor->GetProperty()->SetLineWidth(5);
             medialsheetcrestcurveactor->GetProperty()->SetColor(1,1,0);
             renderer->AddActor(medialsheetcrestcurveactor);

             vtkSmartPointer< vtkSRepVisuPrimitives > visuprimitives = vtkSmartPointer< vtkSRepVisuPrimitives >::New();
             visuprimitives->SetInput(srepfig);
             visuprimitives->Update();

             vtkActorCollection* actorcollection = visuprimitives->GetOuputActors();

             for(unsigned i = 0; i < 4; i++){

                 vtkActor* actor = (vtkActor*) actorcollection->GetItemAsObject(i);

                 renderer->AddActor(actor);
             }

             vtkActor* actor = (vtkActor*) actorcollection->GetItemAsObject(0);//hubs
             //renderer->AddActor(actor);

             actor = (vtkActor*) actorcollection->GetItemAsObject(4);//normals to the sheet
             renderer->AddActor(actor);

             actor = (vtkActor*) actorcollection->GetItemAsObject(5);//uderivative
             //renderer->AddActor(actor);

             actor = (vtkActor*) actorcollection->GetItemAsObject(6);//vderivative
             //renderer->AddActor(actor);



             vtkSRep::VectorIdsType crestids = srepfig->GetCrestMedialAtomsIds();

             //renderWindowInteractor->Start();

             /*//int n = 0;
             //while(n < 1){
                 //n++;
                 for(unsigned i = 0; i < 5; i++){


                     //vtkSmartPointer<vtkSRepInterpolateCrestSpokes> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokes>::New();
                     vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
                     interpolatecrestspokes->SetInput(srepfig);
                     interpolatecrestspokes->SetInterpolationLevel(interpolationlevel);
                     interpolatecrestspokes->SetAtomId(4);
                     if(i == 0){
                         interpolatecrestspokes->SetGamma_t(1);
                         interpolatecrestspokes->SetGamma_theta(1);
                     }else if(i == 1){
                         interpolatecrestspokes->SetGamma_t(0.5);
                         interpolatecrestspokes->SetGamma_theta(0.5);
                     }else if(i == 2){
                         interpolatecrestspokes->SetGamma_t(0.5);
                         interpolatecrestspokes->SetGamma_theta(0.5);
                         interpolatecrestspokes->SetSpokeType(vtkSRep::TOP_SPOKE);
                     }else if(i == 3){
                         interpolatecrestspokes->SetGamma_t(0.5);
                         interpolatecrestspokes->SetGamma_theta(0.5);
                         interpolatecrestspokes->SetSpokeType(vtkSRep::BOTTOM_SPOKE);
                     }else if(i == 4){
                         interpolatecrestspokes->SetGamma_t(0.5);
                         interpolatecrestspokes->SetGamma_theta(0.5);
                         interpolatecrestspokes->SetSpokeType(vtkSRep::CREST_SPOKE);
                     }

                     interpolatecrestspokes->Update();

                     vtkSmartPointer<vtkSRep> srepcrest = interpolatecrestspokes->GetSRepOutput();

                     vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();
                     vtkSmartPointer<vtkCellArray> cellarraycrestspokes = vtkSmartPointer<vtkCellArray>::New();
                     vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();

                     for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){

                         double point[3];

                         srepcrest->GetPoint(i, point);

                         vtkIdType id0 = pointscrestspokes->InsertNextPoint(point[0], point[1], point[2]);

                         vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);
                         vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);

                         for(unsigned j = 0; j < currentspokes.size(); j++){

                             vtkSRep::VNLType p1 = currentspokes[j]*radius[j];

                             vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
                             crestspokeline->GetPointIds()->SetId(0, id0);
                             crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(point[0] + p1[0], point[1] + p1[1], point[2] + p1[2]));

                             cellarraycrestspokes->InsertNextCell(crestspokeline);

                         }

                     }

                     polycrestspokes->SetPoints(pointscrestspokes);
                     polycrestspokes->SetLines(cellarraycrestspokes);

                     vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                     crestspokesmapper->SetInput(polycrestspokes);
                     vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
                     spokesactor->SetMapper(crestspokesmapper);
                     spokesactor->GetProperty()->SetLineWidth(5);
                     spokesactor->GetProperty()->SetColor(0,1,0);
                     renderer->AddActor(spokesactor);
                     renderWindow->Render();


                     renderWindowInteractor->Start();

                     renderer->RemoveActor(spokesactor);

                 }
             //}*/


             //vtkSmartPointer<vtkSRepInterpolateCrestSpokes> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokes>::New();
             vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic> interpolatecrestspokes = vtkSmartPointer<vtkSRepInterpolateCrestSpokesQuartic>::New();
             interpolatecrestspokes->SetInput(srepfig);
             interpolatecrestspokes->SetInterpolationLevel(interpolationlevel);
             interpolatecrestspokes->SetAtomId(9);
             //interpolatecrestspokes->SetGamma_t(0.5);
             //interpolatecrestspokes->SetGamma_theta(0.5);
             //interpolatecrestspokes->SetSpokeType(vtkSRep::CREST_SPOKE);
             interpolatecrestspokes->Update();

             vtkSmartPointer<vtkSRep> srepcrest = interpolatecrestspokes->GetSRepOutput();

             vtkSmartPointer<vtkPolyData> polycrestspokes = vtkSmartPointer<vtkPolyData>::New();
             vtkSmartPointer<vtkCellArray> cellarraycrestspokes = vtkSmartPointer<vtkCellArray>::New();
             vtkSmartPointer<vtkPoints> pointscrestspokes  = vtkSmartPointer<vtkPoints>::New();

             for(unsigned i = 0; i < srepcrest->GetNumberOfPoints(); i++){

                 double point[3];

                 srepcrest->GetPoint(i, point);

                 vtkIdType id0 = pointscrestspokes->InsertNextPoint(point[0], point[1], point[2]);

                 vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);
                 vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);

                 for(unsigned j = 0; j < currentspokes.size(); j++){

                     vtkSRep::VNLType p1 = currentspokes[j]*radius[j];

                     vtkSmartPointer<vtkLine> crestspokeline = vtkSmartPointer<vtkLine>::New();
                     crestspokeline->GetPointIds()->SetId(0, id0);
                     crestspokeline->GetPointIds()->SetId(1, pointscrestspokes->InsertNextPoint(point[0] + p1[0], point[1] + p1[1], point[2] + p1[2]));

                     cellarraycrestspokes->InsertNextCell(crestspokeline);

                 }

             }

             polycrestspokes->SetPoints(pointscrestspokes);
             polycrestspokes->SetLines(cellarraycrestspokes);

             vtkSmartPointer<vtkPolyDataMapper> crestspokesmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
             crestspokesmapper->SetInput(polycrestspokes);
             vtkSmartPointer<vtkActor>  spokesactor = vtkActor::New();
             spokesactor->SetMapper(crestspokesmapper);
             spokesactor->GetProperty()->SetLineWidth(5);
             spokesactor->GetProperty()->SetColor(0,1,0);
             renderer->AddActor(spokesactor);


             vtkPolyData* interpolatedcrest = interpolatecrestspokes->GetOutput();

              vtkSmartPointer<vtkPolyDataMapper> crestspokescurvemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
              crestspokescurvemapper->SetInputConnection(interpolatedcrest->GetProducerPort());
              vtkSmartPointer<vtkActor>  crestspokesactor = vtkActor::New();
              crestspokesactor->SetMapper(crestspokescurvemapper);
              crestspokesactor->GetProperty()->SetLineWidth(5);
              crestspokesactor->GetProperty()->SetColor(1,0,0);
              renderer->AddActor(crestspokesactor);




              /*vtkSmartPointer<vtkAppendPolyData> appendpoly = vtkSmartPointer<vtkAppendPolyData>::New();
              appendpoly->AddInput(interpolatedcrest);
              appendpoly->AddInput(interpolatedmedialspokes);
              appendpoly->Update();

              vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
              clean->SetInput(appendpoly->GetOutput());
              clean->Update();


              vtkSmartPointer<vtkPolyDataMapper> srepmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
              srepmapper->SetInput(clean->GetOutput());
              vtkSmartPointer<vtkActor>  srepactor = vtkActor::New();
              srepactor->SetMapper(srepmapper);
              srepactor->GetProperty()->SetColor(0,0.5,0.5);
              renderer->AddActor(srepactor);






              vtkSmartPointer<vtkSTLWriter> polywriter = vtkSmartPointer<vtkSTLWriter>::New();
              polywriter->SetInput(appendpoly->GetOutput());
              polywriter->SetFileName("srepout.vtk");
              polywriter->Write();*/


             renderWindow->Render();
             renderWindowInteractor->Start();
             renderer->RemoveAllViewProps();
         }



    return 0;
}








