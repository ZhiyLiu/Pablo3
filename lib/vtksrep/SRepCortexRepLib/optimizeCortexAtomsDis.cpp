#include "optimizeCortexAtomsDis.h"

#include "iostream"

#include "vtkIdList.h"
#include "vtkCell.h"

using namespace std;

optimizeCortexAtomsDis::optimizeCortexAtomsDis()
    : vnl_least_squares_function(3,4,no_gradient)
{
    m_NPoints = 9;
}


optimizeCortexAtomsDis::~optimizeCortexAtomsDis()
{
}

void optimizeCortexAtomsDis::f(vnl_vector< double > const &x, vnl_vector< double > &fx){

    //cout<<x<<endl;
    //fx.fill(0);

    //if(fabs(100.0 - x.magnitude()) > 0.25){
    //    fx.fill(10000);
    //}else{

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

                    fx[i] = pow((x - vnlnextpoint).magnitude(), 2)*pow(10.0,0);
                }
            }
        }
    //}
}
