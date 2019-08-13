#ifndef vtkOptimizeCortexAtomsCurv_H
#define vtkOptimizeCortexAtomsCurv_H

#include "vtkAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include <optimizeCortexAtoms.h>

#include "vtksrep.h"

#include <vector>

using namespace std;

class vtkOptimizeCortexAtomsCurv : public vtkAlgorithm
{
public:
    static vtkOptimizeCortexAtomsCurv *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< VectorVNLType > VectorVectorVNLType;


    void SetInputCortexSurf(vtkPolyData* pial, vtkPolyData* orig){
        m_Pial = pial;
        m_Orig = orig;
    }

    void SetInputSphere(vtkPolyData* sphere){
        m_Sphere = sphere;
    }

    void SetInput(vtkSRep* srep){
        m_SRep = srep;
    }

    void SetInputCurvatureArray(vtkDataArray* curvatureArray){
        m_CurvatureArray = curvatureArray;
    }

    vtkSRep* GetOutput(){
        return m_SRep;
    }




    void Update();

protected:
    vtkOptimizeCortexAtomsCurv();

    ~vtkOptimizeCortexAtomsCurv();


private:


    vtkSRep* m_SRep;
    vtkPolyData* m_Sphere;
    vtkPolyData* m_Pial;
    vtkPolyData* m_Orig;
    int m_Niter;
    vtkDataArray* m_CurvatureArray;

};

#endif
