#include "optimizeCortexAtoms.h"

#include "iostream"

#include "vtkIdList.h"
#include "vtkCell.h"

using namespace std;

OptimizeCortexAtoms::OptimizeCortexAtoms()
    : vnl_least_squares_function(3,13,no_gradient)
{
    m_Locator = 0;
}


OptimizeCortexAtoms::~OptimizeCortexAtoms()
{
}

void OptimizeCortexAtoms::f(vnl_vector< double > const &x, vnl_vector< double > &fx){

    //cout<<x<<endl;
    double point[3];
    point[0] = x[0];
    point[1] = x[1];
    point[2] = x[2];


    vtkSmartPointer<vtkIdList> pointsid = vtkSmartPointer<vtkIdList>::New();
    m_Locator->FindClosestNPoints(9, point, pointsid);

    /*vnl_vector<double> avgpoint(x);
    avgpoint.fill(0);
    for(unsigned i = 0; i < pointsid->GetNumberOfIds(); i++){

        double current[3];
        m_Sphere->GetPoint(pointsid->GetId(i), current);
        vnl_vector<double> temp(current, 3);

        avgpoint += temp;

    }
    avgpoint/=pointsid->GetNumberOfIds();*/

    //double totalcurv = 0;

    //if(sqrt(pow(m_X0[0] - avgpoint[0], 2) + pow(m_X0[1] - avgpoint[1], 2) + pow(m_X0[2] - avgpoint[2], 2)) >= m_MaxMov){
    if(x.magnitude() > 110 || x.magnitude() < 90 || (m_X0 - x).magnitude() > m_MaxMov){
        fx.fill(VTK_DOUBLE_MAX);
    }else{

        double avgcurv = 0;
        double avgweight = 0;
        double ro = -1;

        for(unsigned i = 0; i < pointsid->GetNumberOfIds(); i++){

            double curvtuple[1];
            m_Curvature->GetTuple(pointsid->GetId(i), curvtuple);

            avgcurv += curvtuple[0];

            double current[3];
            m_Sphere->GetPoint(pointsid->GetId(i), current);

            double dist = sqrt(pow(current[0] - x[0], 2) + pow(current[1] - x[1], 2) + pow(current[2] - x[2], 2));
            avgweight += pow(dist,ro);

        }

        avgcurv/=pointsid->GetNumberOfIds();
        avgweight/=pointsid->GetNumberOfIds();

        for(unsigned i = 0; i < pointsid->GetNumberOfIds(); i++){

            double current[3];
            m_Sphere->GetPoint(pointsid->GetId(i), current);

            double dist = sqrt(pow(current[0] - x[0], 2) + pow(current[1] - x[1], 2) + pow(current[2] - x[2], 2));
            double weight = pow(dist,ro);

            double curvtuple[1];
            m_Curvature->GetTuple(pointsid->GetId(i), curvtuple);
            //double val = (curvtuple[0]-avgcurv);
            //double val = weight*(1.0 - fabs(curvtuple[0]))/avgweight;
            double val = weight*(fabs(curvtuple[0]))/avgweight;

            fx[i] = val;
        }
        //cout<<fx<<endl;
    }

    vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
    m_SRep->GetPointCells(m_Id, idlist);

    for(unsigned i = 0; i < idlist->GetNumberOfIds(); i++){
        vtkCell* cell = m_SRep->GetCell(idlist->GetId(i));
        for(unsigned j = 0; j < cell->GetPointIds()->GetNumberOfIds(); j++){

            if(m_Id == cell->GetPointId(j)){
                double nextpoint[3];

                if( j < 3){
                    m_SRep->GetPoint(cell->GetPointId(j+1), nextpoint);
                }else{
                    m_SRep->GetPoint(cell->GetPointId(0), nextpoint);
                }

                vnl_vector<double> vnlnextpoint(nextpoint, 3);

                fx[9+i] = pow((x - vnlnextpoint).magnitude(), 2)*pow(10,-0.5);
            }
        }
    }
}
