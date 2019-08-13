#include "vtkoptimizeCortexAtomsCurv.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkLine.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"


#include "optimizeCortexAtomsDis.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkOptimizeCortexAtomsCurv);



vtkOptimizeCortexAtomsCurv::vtkOptimizeCortexAtomsCurv()
{
    m_Sphere = 0;
    m_Pial = 0;
    m_Orig = 0;
    m_Niter = 1;

}

vtkOptimizeCortexAtomsCurv::~vtkOptimizeCortexAtomsCurv()
{
}


// Superclass method to update the pipeline
void vtkOptimizeCortexAtomsCurv::Update(){

    vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
    locator->SetDataSet(m_Sphere);
    locator->AutomaticOn();
    locator->SetNumberOfPointsPerBucket(3);
    locator->BuildLocator();

    int niter = 0;
    while(niter < m_Niter){
        niter++;
        for(unsigned i = 0; i < m_SRep->GetNumberOfPoints(); i++){

            double point[3];
            m_SRep->GetPoint(i, point);

            vnl_vector<double> x0(point, 3);

            vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
            m_SRep->GetPointCells(i, cellids);
            double maxmov = VTK_DOUBLE_MAX;

            for(unsigned j = 0; j < cellids->GetNumberOfIds(); j++){
                vtkCell* quad = m_SRep->GetCell(cellids->GetId(j));
                for(unsigned k = 0; k < quad->GetNumberOfPoints(); k++){
                    if(quad->GetPointIds()->GetId(k) != i){
                        vnl_vector<double> temp(quad->GetPoints()->GetPoint(k), 3);
                        double tempmaxmov = (x0 - temp).magnitude();
                        if(tempmaxmov < maxmov){
                            maxmov = tempmaxmov;
                        }
                    }
                }
            }
            //cout<<maxmov<<endl;

            OptimizeCortexAtoms optimizecortex;

            optimizecortex.SetCurvatureArray(m_CurvatureArray);
            optimizecortex.SetSphere(m_Sphere);
            optimizecortex.SetX0(x0);
            optimizecortex.SetPointLocator(locator);
            optimizecortex.SetMaxMov(maxmov);
            optimizecortex.SetSRep(m_SRep);
            optimizecortex.SetId(i);

            vnl_levenberg_marquardt levenberg(optimizecortex);

            //cout<<endl<<"id = "<<i<<endl;

            //levenberg.set_max_function_evals(20);
            levenberg.minimize(x0);
            //cout<<endl;
            vtkSmartPointer<vtkIdList> pointsid = vtkSmartPointer<vtkIdList>::New();
            locator->FindClosestNPoints(3, x0.data_block(), pointsid);

            vnl_vector<double> avgpointpial(3);
            avgpointpial.fill(0);
            vnl_vector<double> avgpointorig(3);
            avgpointorig.fill(0);
            for(unsigned j = 0; j < pointsid->GetNumberOfIds(); j++){
                vtkIdType pid = pointsid->GetId(j);
                double closestpoint[3];

                if(niter == 1 && m_Orig){
                    m_Orig->GetPoint(pid, closestpoint);
                }else{
                    m_Sphere->GetPoint(pid, closestpoint);
                }
                vnl_vector<double> vnlorig(closestpoint, 3);
                avgpointorig += vnlorig;

                if(niter == m_Niter && m_Pial){
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
