#include "vtkoptimizeCortexSpokes.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkVertex.h"
#include "vtkCleanPolyData.h"
#include "vtkSmoothPolyDataFilter.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"

#include "vtkinterpolatecurve.h"

#include "optimizeCortexSpokes.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkOptimizeCortexSpokes);

#include <map>

vtkOptimizeCortexSpokes::vtkOptimizeCortexSpokes()
{
    m_Sphere = 0;
    m_Pial = 0;
    m_Orig = 0;

}

vtkOptimizeCortexSpokes::~vtkOptimizeCortexSpokes()
{
}


// Superclass method to update the pipeline
void vtkOptimizeCortexSpokes::Update(){

    if(m_SamplePoints.size() > 1){


        vtkSmartPointer< vtkInterpolateCurve > interpolatecurve = vtkSmartPointer< vtkInterpolateCurve >::New();
        interpolatecurve->SetCurvePoints(m_SamplePoints);
        interpolatecurve->SetCyclic(false);
        interpolatecurve->SetInterpolationLevel(3);
        interpolatecurve->Update();

        vtkSRep::VectorIdsType srependatomsid = m_SRep->GetCrestMedialAtomsIds();


        double step = 1.0/(srependatomsid.size() - 1);
        for(double t = 0, n = 0; t < 1; t+=step, n++){
            vtkInterpolateCurve::VNLType point = interpolatecurve->EvaluateFunction(t);
            m_SRep->GetPoints()->SetPoint(srependatomsid[n], point[0], point[1], point[2]);
        }


        /*vector< vnl_vector< double > > samplepoints;
        samplepoints.insert(samplepoints.begin(), m_SamplePoints.begin() + m_SamplePoints.size()/2.0, m_SamplePoints.end());

        vtkSmartPointer< vtkInterpolateCurve > interpolatecurve = vtkSmartPointer< vtkInterpolateCurve >::New();
        interpolatecurve->SetCurvePoints(samplepoints);
        interpolatecurve->SetCyclic(false);
        interpolatecurve->SetInterpolationLevel(3);
        interpolatecurve->Update();

        vtkSRep::VectorIdsType srependatomsid = m_SRep->GetCrestMedialAtomsIds(3);


        double step = 1.0/(srependatomsid.size() - 1);
        for(double t = 0, n = 0; t < 1; t+=step, n++){
            vtkInterpolateCurve::VNLType point = interpolatecurve->EvaluateFunction(t);
            m_SRep->GetPoints()->SetPoint(srependatomsid[n], point[0], point[1], point[2]);
        }

        samplepoints.clear();
        samplepoints.insert(samplepoints.begin(), m_SamplePoints.begin(), m_SamplePoints.begin() + m_SamplePoints.size()/2.0);
        interpolatecurve = vtkSmartPointer< vtkInterpolateCurve >::New();
        interpolatecurve->SetCurvePoints(samplepoints);
        interpolatecurve->SetCyclic(false);
        interpolatecurve->SetInterpolationLevel(3);
        interpolatecurve->Update();


        srependatomsid.clear();
        vtkSRep::VectorIdsType srependatomsid0 = m_SRep->GetCrestMedialAtomsIds(0);
        vtkSRep::VectorIdsType srependatomsid1 = m_SRep->GetCrestMedialAtomsIds(1);
        vtkSRep::VectorIdsType srependatomsid2 = m_SRep->GetCrestMedialAtomsIds(2);
        srependatomsid.insert(srependatomsid.end(), srependatomsid0.begin(), srependatomsid0.end());
        srependatomsid.insert(srependatomsid.end(), srependatomsid1.begin(), srependatomsid1.end());
        srependatomsid.insert(srependatomsid.end(), srependatomsid2.begin(), srependatomsid2.end());


        step = 1.0/(srependatomsid.size() - 1);
        for(double t = 0, n = 0; t < 1; t+=step, n++){
            vtkInterpolateCurve::VNLType point = interpolatecurve->EvaluateFunction(t);
            m_SRep->GetPoints()->SetPoint(srependatomsid[n], point[0], point[1], point[2]);
        }*/

    }




    if(m_Sphere){
        vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
        locator->SetDataSet(m_Sphere);
        locator->AutomaticOn();
        locator->SetNumberOfPointsPerBucket(3);
        locator->BuildLocator();

        for(unsigned i = 0; i < m_SRep->GetNumberOfPoints(); i++){

            double point[3];
            m_SRep->GetPoint(i, point);

            vtkSmartPointer<vtkIdList> pointsid = vtkSmartPointer<vtkIdList>::New();
            locator->FindClosestNPoints(3, point, pointsid);

            vnl_vector<double> avgpointpial(3);
            avgpointpial.fill(0);
            vnl_vector<double> avgpointorig(3);
            avgpointorig.fill(0);
            for(unsigned j = 0; j < pointsid->GetNumberOfIds(); j++){
                vtkIdType pid = pointsid->GetId(j);
                double closestpoint[3];

                if(m_Orig){
                    m_Orig->GetPoint(pid, closestpoint);
                }else{
                    m_Sphere->GetPoint(pid, closestpoint);
                }
                vnl_vector<double> vnlorig(closestpoint, 3);
                avgpointorig += vnlorig;

                if(m_Pial){
                    m_Pial->GetPoint(pid, closestpoint);
                }else{
                    m_Sphere->GetPoint(pid, closestpoint);
                }
                vnl_vector<double> vnlpial(closestpoint, 3);
                avgpointpial += vnlpial;
            }

            avgpointpial /= pointsid->GetNumberOfIds();
            avgpointorig /= pointsid->GetNumberOfIds();

            vnl_vector<double> midpoint = (avgpointpial + avgpointorig)/2.0;

            vtkSRep::VNLType s(3);
            vtkSRep::VNLType s1(3);

            s = avgpointpial - midpoint;
            s1 = midpoint - avgpointorig;

            m_SRep->SetSpoke(i, vtkSRep::TOP_SPOKE, s);
            m_SRep->SetSpoke(i, vtkSRep::BOTTOM_SPOKE, s1);
            m_SRep->GetPoints()->SetPoint(i, midpoint[0], midpoint[1], midpoint[2]);
        }
    }
}
