#ifndef vtkOptimizeSRepSampling_H
#define vtkOptimizeSRepSampling_H

#include "vtkAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "optimizeCortexAtomsDis.h"

#include "vtksrep.h"

#include <vector>

using namespace std;

class vtkOptimizeSRepSampling : public vtkAlgorithm
{
public:
    static vtkOptimizeSRepSampling *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< VectorVNLType > VectorVectorVNLType;

    void SetInput(vtkSRep* srep){
        m_SRep = srep;
    }

    vtkSRep* GetOutput(){
        return m_SRep;
    }

    virtual void Update();

    void SetMax(double max){
        m_Max = max;
    }
    void SetMin(double min){
        m_Min = min;
    }
    void SetNumP(int nump){
        m_NumP = nump;
    }

    void SetSphere(vtkPolyData* sphere){
        m_Sphere = sphere;
    }

    void SetRegisteredPoints(VectorVNLType registered){
        m_RegisteredPoints = registered;
    }

protected:
    vtkOptimizeSRepSampling();

    ~vtkOptimizeSRepSampling();


private:


    vtkSmartPointer<vtkSRep> m_SRep;
    double m_Max;
    double m_Min;
    int m_NumP;
    vtkPolyData* m_Sphere;
    double m_SphereRadius;
    VectorVNLType m_RegisteredPoints;



};

#endif
