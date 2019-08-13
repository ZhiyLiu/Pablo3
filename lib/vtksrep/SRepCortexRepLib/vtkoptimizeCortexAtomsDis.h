#ifndef vtkOptimizeCortexAtomsDis_H
#define vtkOptimizeCortexAtomsDis_H

#include "vtkAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "optimizeCortexAtomsDis.h"

#include "vtksrep.h"

#include <vector>

using namespace std;

class vtkOptimizeCortexAtomsDis : public vtkAlgorithm
{
public:
    static vtkOptimizeCortexAtomsDis *New();


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

    vtkSRep* GetOutput(){
        return m_SRep;
    }

    void SetKeepSpokes(bool keep){
        m_KeepSpokes = keep;
    }

    void SetLabelArray(vtkDataArray* labelarray){
        m_LabelArray = labelarray;
    }

    void Update();

protected:
    vtkOptimizeCortexAtomsDis();

    ~vtkOptimizeCortexAtomsDis();


private:


    vtkSRep* m_SRep;
    vtkPolyData* m_Sphere;
    vtkPolyData* m_Pial;
    vtkPolyData* m_Orig;
    vtkDataArray* m_LabelArray;
    bool m_KeepSpokes;

};

#endif
