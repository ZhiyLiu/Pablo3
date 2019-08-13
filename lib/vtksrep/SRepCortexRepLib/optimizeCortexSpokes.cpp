#include "optimizeCortexSpokes.h"

#include "iostream"

#include "vtkIdList.h"



using namespace std;

optimizeCortexSpokes::optimizeCortexSpokes()
    : vnl_least_squares_function(3,64,no_gradient)
{

}


optimizeCortexSpokes::~optimizeCortexSpokes()
{
}

void optimizeCortexSpokes::f(vnl_vector< double > const &x, vnl_vector< double > &fx){

    //cout<<x<<endl;    
    vnl_vector<double> s0(3);
    s0[0] = x[0];
    s0[1] = x[1];
    s0[2] = x[2];
    s0.normalize();

    /*vnl_vector<double> s1(3);
    s1[0] = x[3];
    s1[1] = x[4];
    s1[2] = x[5];
    s1.normalize();*/

    cout<<s0<<endl;
    //cout<<s1<<endl;

    double r0 = m_Srep->GetSpokeRadius(m_AtomId, vtkSRep::TOP_SPOKE);
    //double r1 = m_Srep->GetSpokeRadius(m_AtomId, vtkSRep::BOTTOM_SPOKE);
    m_Srep->SetSpoke(m_AtomId, vtkSRep::TOP_SPOKE, s0*r0);
    //m_Srep->SetSpoke(m_AtomId, vtkSRep::BOTTOM_SPOKE, s1*r1);
    fx.fill(0);

    {

        m_Medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();
        m_Medialspokesinterpolator->SetInput(m_Srep);
        m_Medialspokesinterpolator->SetInterpolationLevel(2);
        m_Medialspokesinterpolator->SetAtomId(m_AtomId);
        m_Medialspokesinterpolator->SetSpokeType(vtkSRep::TOP_SPOKE);
        m_Medialspokesinterpolator->Update();
        vtkSmartPointer<vtkPolyData> interpolatedsurf = m_Medialspokesinterpolator->GetOutput();

        //cout<<interpolatedsurf->GetNumberOfPoints()<<endl;
        for(unsigned i = 0; i < interpolatedsurf->GetNumberOfPoints() && i < fx.size(); i++){
            double point[3];
            interpolatedsurf->GetPoint(i, point);

            vtkIdType id = m_Locator0->FindClosestPoint(point);
            double closestp[3];
            m_Locator0->GetDataSet()->GetPoint(id, closestp);

            vnl_vector<double> p0(point, 3);
            vnl_vector<double> p1(closestp, 3);
            fx[i] = (p0 - p1).magnitude();
        }

    }
    /*if(false){
        vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > medialspokesinterpolator = vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite >::New();
        medialspokesinterpolator->SetInput(m_Srep);
        medialspokesinterpolator->SetInterpolationLevel(3);
        medialspokesinterpolator->SetAtomId(m_AtomId);
        medialspokesinterpolator->SetSpokeType(vtkSRep::BOTTOM_SPOKE);
        medialspokesinterpolator->Update();
        vtkPolyData* interpolatedsurf = medialspokesinterpolator->GetOutput();

        for(unsigned i = 0; i < interpolatedsurf->GetNumberOfPoints() && (i+256) < fx.size(); i++){
            double point[3];
            interpolatedsurf->GetPoint(i, point);

            vtkIdType id = m_Locator1->FindClosestPoint(point);
            double closestp[3];
            m_Locator1->GetDataSet()->GetPoint(id, closestp);

            vnl_vector<double> p0(point, 3);
            vnl_vector<double> p1(closestp, 3);
            fx[i+64] = (p0 - p1).magnitude();
        }
    }*/

    //cout<<fx.magnitude()<<endl;


}
