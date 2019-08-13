#ifndef vtkOptimizeCortexSpokes_H
#define vtkOptimizeCortexSpokes_H

#include "vtkAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "optimizeCortexAtomsDis.h"

#include "vtksrep.h"

#include <vector>

using namespace std;

class vtkOptimizeCortexSpokes : public vtkAlgorithm
{
public:
    static vtkOptimizeCortexSpokes *New();


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

    void SetInputLabelStructure(int label){
        m_LabelStructure = label;
    }

    void SetInputVectorLabel(vnl_vector<int> vectlabel){
        m_VectLabel = vectlabel;
    }

    void Update();

    vtkPolyData* GetOuputCorpusContour(){
        return m_CorpusPoly;
    }

    void SetInputSamplePoints(vector< vnl_vector< double > > samplepoints){
        m_SamplePoints = samplepoints;
    }

protected:
    vtkOptimizeCortexSpokes();

    ~vtkOptimizeCortexSpokes();


private:


    vtkSRep* m_SRep;
    vtkPolyData* m_Sphere;
    vtkPolyData* m_Pial;
    vtkPolyData* m_Orig;
    int m_LabelStructure;
    vnl_vector<int> m_VectLabel;
    vtkSmartPointer< vtkPolyData > m_CorpusPoly;
    vector< vnl_vector< double > > m_SamplePoints;

};

#endif
