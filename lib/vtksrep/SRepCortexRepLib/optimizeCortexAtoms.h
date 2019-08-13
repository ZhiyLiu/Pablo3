#ifndef OPTIMIZECORTEXATOMS_H
#define OPTIMIZECORTEXATOMS_H

#include "vnl/vnl_least_squares_function.h"

#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkPointLocator.h>
#include <vtkSmartPointer.h>
#include <vtksrep.h>

class OptimizeCortexAtoms : public vnl_least_squares_function
{
public:
    OptimizeCortexAtoms();    
    ~OptimizeCortexAtoms();

    virtual void f(vnl_vector< double > const &x, vnl_vector< double > &fx);


    vnl_vector< double > GetCoefficients(double q0, double thetamax = 1);


    void SetX0(vnl_vector<double> x0){
        m_X0 = x0;
    }
    void SetMaxMov(double maxmov){
        m_MaxMov = maxmov;
    }

    void SetPointLocator(vtkSmartPointer<vtkPointLocator> locator){
        m_Locator = locator;
    }

    void SetSphere(vtkSmartPointer<vtkPolyData> sphere){
        m_Sphere = sphere;
    }

    void SetCurvatureArray(vtkSmartPointer<vtkDataArray> curvature){
        m_Curvature = curvature;
    }

    void SetId(vtkIdType id){
        m_Id = id;
    }
    void SetSRep(vtkSRep* srep){
        m_SRep = srep;
    }

private:

    vnl_vector<double> m_X0;
    double m_MaxMov;
    vtkSmartPointer<vtkPointLocator> m_Locator;
    vtkSmartPointer<vtkPolyData> m_Sphere;
    vtkSmartPointer<vtkDataArray> m_Curvature;
    vtkIdType m_Id;
    vtkSRep* m_SRep;

};

#endif // OPTIMIZECORTEXATOMS_H
