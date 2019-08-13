#include "vtkoptimizeCortexAtomsDis.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkLine.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"


#include "optimizeCortexAtomsDis.h"
#include "optimizeCortexSpokes.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkOptimizeCortexAtomsDis);

#include "vtkPolyDataNormals.h"

#include "vtkPointData.h"
#include "map"

vtkOptimizeCortexAtomsDis::vtkOptimizeCortexAtomsDis()
{
    m_Sphere = 0;
    m_Pial = 0;
    m_Orig = 0;
    m_KeepSpokes = false;
    m_LabelArray = 0;

}

vtkOptimizeCortexAtomsDis::~vtkOptimizeCortexAtomsDis()
{
}


// Superclass method to update the pipeline
void vtkOptimizeCortexAtomsDis::Update(){


    /*vtkSmartPointer<vtkPointLocator> locatorsrep = vtkSmartPointer<vtkPointLocator>::New();
    locatorsrep->SetDataSet(m_SRep);
    locatorsrep->AutomaticOn();
    locatorsrep->SetNumberOfPointsPerBucket(3);
    locatorsrep->BuildLocator();

    for(unsigned i = 0; i < m_SRep->GetNumberOfPoints();i++){

        vtkSmartPointer<vtkIdList> cellsid = vtkSmartPointer<vtkIdList>::New();
        m_SRep->GetPointCells(i, cellsid);

        if(cellsid->GetNumberOfIds() == 4){
            double point[3];
            m_SRep->GetPoint(i, point);

            vtkSRep::VectorVNLType spokes = m_SRep->GetSpokes(i);
            if(m_KeepSpokes){
                vtkSRep::VNLType vnlpoint(point, 3);
                for(unsigned j = 0; j < spokes.size(); j++){
                    spokes[j] += vnlpoint;
                    //cout<<spokes[j]<<endl;
                }
            }

            optimizeCortexAtomsDis optimizecortex;
            optimizecortex.SetPointLocator(locatorsrep);
            optimizecortex.SetSRep(m_SRep);
            optimizecortex.SetId(i);

            vnl_levenberg_marquardt levenberg(optimizecortex);

            vnl_vector<double> x0(point, 3);

            levenberg.minimize(x0);

            if(m_KeepSpokes){
                for(unsigned j = 0; j < spokes.size(); j++){
                    spokes[j] -= x0;
                    cout<<spokes[j]<<endl;
                }
                m_SRep->SetSpokes(i, spokes);
            }

            m_SRep->GetPoints()->SetPoint(i, x0[0], x0[1], x0[2]);

        }
    }

    for(unsigned i = 0; i < m_SRep->GetNumberOfPoints();i++){

        vtkSmartPointer<vtkIdList> cellsid = vtkSmartPointer<vtkIdList>::New();
        m_SRep->GetPointCells(i, cellsid);

        if(cellsid->GetNumberOfIds() < 4){
            double point[3];
            m_SRep->GetPoint(i, point);

            optimizeCortexAtomsDis optimizecortex;
            optimizecortex.SetPointLocator(locatorsrep);
            optimizecortex.SetSRep(m_SRep);
            optimizecortex.SetId(i);

            vnl_levenberg_marquardt levenberg(optimizecortex);

            vnl_vector<double> x0(point, 3);

            levenberg.minimize(x0);


            m_SRep->GetPoints()->SetPoint(i, x0[0], x0[1], x0[2]);

        }
    }*/

    if(m_Sphere){

        vtkSmartPointer< vtkPolyData > midsurf = 0;

        if(m_Pial && m_Orig){

            midsurf = vtkSmartPointer< vtkPolyData >::New();
            vtkSmartPointer< vtkPoints > midsurfpoints = vtkSmartPointer< vtkPoints >::New();
            for(unsigned i = 0; i < m_Pial->GetNumberOfPoints(); i++){

                double p0[3];
                m_Pial->GetPoint(i, p0);

                double p1[3];
                m_Orig->GetPoint(i, p1);

                double p2[3];
                p2[0] = (p1[0] + p2[0])/2.0;
                p2[1] = (p1[1] + p2[1])/2.0;
                p2[2] = (p1[2] + p2[2])/2.0;

                midsurfpoints->InsertNextPoint(p2[0], p2[1], p2[2]);

            }
            midsurf->SetPoints(midsurfpoints);
            midsurf->SetPolys(m_Pial->GetPolys());

            vtkSmartPointer<vtkPolyDataNormals> polynormals = vtkSmartPointer<vtkPolyDataNormals>::New();
            polynormals->SetInput(midsurf);
            polynormals->ComputePointNormalsOn();
            polynormals->Update();
            midsurf = polynormals->GetOutput();
        }


        vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
        locator->SetDataSet(m_Sphere);
        locator->AutomaticOn();
        locator->SetNumberOfPointsPerBucket(3);
        locator->BuildLocator();

        vtkSmartPointer< vtkDataArray > labelarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
        vtkSmartPointer< vtkDataArray > normalarray = (vtkDataArray*)vtkDataArray::CreateArray(VTK_DOUBLE);
        normalarray->SetNumberOfComponents(3);

        for(unsigned i = 0; i < m_SRep->GetNumberOfPoints(); i++){

            double point[3];
            m_SRep->GetPoint(i, point);

            vtkSmartPointer<vtkIdList> pointsid = vtkSmartPointer<vtkIdList>::New();
            locator->FindClosestNPoints(3, point, pointsid);

            vnl_vector<double> vnlpoint(point, 3);
            vnl_vector<double> avgpointpial(3);
            avgpointpial.fill(0);
            vnl_vector<double> avgpointorig(3);
            avgpointorig.fill(0);

            vnl_vector<double> avgnormal(3);
            avgnormal.fill(0);

            double avgweight = 0;
            double ro = -0.8;
            for(unsigned j = 0; j < pointsid->GetNumberOfIds(); j++){
                double pointinf[3];

                vtkIdType pid = pointsid->GetId(j);

                m_Sphere->GetPoint(pid, pointinf);
                vtkSRep::VNLType vnlpointinf(pointinf, 3);
                double dist = (vnlpoint - vnlpointinf).magnitude();
                avgweight += pow(dist, ro);
            }
            //avgweight /= pointsid->GetNumberOfIds();

            double label = 0;
            double maxweight = 0;

            for(unsigned j = 0; j < pointsid->GetNumberOfIds(); j++){
                double pointinf[3];

                vtkIdType pid = pointsid->GetId(j);

                m_Sphere->GetPoint(pid, pointinf);
                vtkSRep::VNLType vnlpointinf(pointinf, 3);
                double dist = (vnlpoint - vnlpointinf).magnitude();
                double weight = pow(dist, ro)/avgweight;

                if(m_Orig){
                    m_Orig->GetPoint(pid, pointinf);
                }else{
                    m_Sphere->GetPoint(pid, pointinf);
                }
                vnl_vector< double > temporig(pointinf, 3);
                avgpointorig += temporig*weight;
                temporig.clear();

                if(m_Pial){
                    m_Pial->GetPoint(pid, pointinf);
                }else{
                    m_Sphere->GetPoint(pid, pointinf);
                }
                vnl_vector< double > temppial(pointinf, 3);
                avgpointpial += temppial*weight;
                temppial.clear();

                if(midsurf){
                    midsurf->GetPointData()->GetNormals()->GetTuple(pid, pointinf);
                }
                vnl_vector< double > tempnormal(pointinf, 3);
                avgnormal += tempnormal*weight;
                tempnormal.clear();

                if(m_LabelArray){
                    if(maxweight < weight){
                        label = m_LabelArray->GetTuple1(pid);
                        maxweight = weight;
                    }
                    //label += m_LabelArray->GetTuple1(pid)*weight;
                }
            }

            labelarray->InsertNextTuple1(label);

            vnl_vector<double> midpoint = (avgpointpial + avgpointorig)/2.0;


            vtkSRep::VNLType s(3);
            vtkSRep::VNLType s1(3);

            s = avgpointpial - midpoint;
            s1 = avgpointorig - midpoint;

            if(midsurf){
                avgnormal.normalize();
                normalarray->InsertNextTuple3(avgnormal[0], avgnormal[1], avgnormal[2]);
            }

            m_SRep->SetSpoke(i, vtkSRep::TOP_SPOKE, s);
            m_SRep->SetSpoke(i, vtkSRep::BOTTOM_SPOKE, s1);
            m_SRep->GetPoints()->SetPoint(i, midpoint[0], midpoint[1], midpoint[2]);
        }


        m_SRep->GetPointData()->SetScalars(labelarray);

        if(midsurf){
            //m_SRep->GetPointData()->SetNormals(normalarray);
        }


        /*vtkSmartPointer<vtkPointLocator> locator0 = vtkSmartPointer<vtkPointLocator>::New();
        locator0->SetDataSet(m_Pial);
        locator0->AutomaticOn();
        locator0->SetNumberOfPointsPerBucket(6);
        locator0->BuildLocator();

        vtkSmartPointer<vtkPointLocator> locator1 = vtkSmartPointer<vtkPointLocator>::New();
        locator1->SetDataSet(m_Orig);
        locator1->AutomaticOn();
        locator1->SetNumberOfPointsPerBucket(6);
        locator1->BuildLocator();



        for(unsigned i = 0; i < m_SRep->GetNumberOfPoints(); i++){


            optimizeCortexSpokes optimizecortex;
            optimizecortex.SetPointLocator0(locator0);
            optimizecortex.SetPointLocator1(locator1);
            optimizecortex.SetSrep(m_SRep);
            optimizecortex.SetAtomId(i);



            vnl_vector<double> x0(6);
            vnl_vector<double> s0 = m_SRep->GetSpoke(i, vtkSRep::TOP_SPOKE, true);
            vnl_vector<double> s1 = m_SRep->GetSpoke(i, vtkSRep::BOTTOM_SPOKE, true);
            x0[0] = s0[0];
            x0[1] = s0[1];
            x0[2] = s0[2];

            x0[3] = s1[0];
            x0[4] = s1[1];
            x0[5] = s1[2];

            vnl_levenberg_marquardt levenberg(optimizecortex);
            levenberg.set_max_function_evals(5);
            levenberg.minimize(s0);

            cout<<i<<endl;

        }*/
    }




}
