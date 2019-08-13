#ifndef vtkSRepInterpolateCrestSpokesQuartic_H
#define vtkSRepInterpolateCrestSpokesQuartic_H

#include "vtkPolyDataAlgorithm.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"
#include "vnl/vnl_numeric_traits.h"

#include "vtksrep.h"
#include "vtkinterpolatecurve.h"
#include "vtkSmartPointer.h"

using namespace std;

class vtkSRepInterpolateCrestSpokesQuartic : public vtkPolyDataAlgorithm
{
public:
    static vtkSRepInterpolateCrestSpokesQuartic *New();


    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< vtkSRep::VectorVNLType > VectorSRepVectorVNLType;
    typedef vtkSRep::VNLMatrixType VNLMatrixType;



    void SetInterpolationLevel(int level){
        m_InterpolationLevel = level;
    }
    int GetInterpolationLevel(){
        return m_InterpolationLevel;
    }

    vtkSmartPointer<vtkSRep> GetSRepOutput(){
        return m_SRepOutput;
    }

    void SetAtomId(vtkIdType atomid){
        m_AtomId = atomid;
    }

    void SetGamma_t(double gammat){
        m_Gamma_t = gammat;
    }

    void SetGamma_theta(double gammatheta){
        m_Gamma_theta = gammatheta;
    }

    void SetSpokeType(vtkIdType spoketype){
        m_SpokeType = spoketype;
    }

    /*
      Set to true for quasi tube interpolation
    */
    void SetUseAllSpokes(bool setall){
        m_UseAllSpokes = setall;
    }
    void SetCyclicSpokes(bool cyclic){
        m_CyclicSpokes = cyclic;
    }

protected:
    vtkSRepInterpolateCrestSpokesQuartic();

    ~vtkSRepInterpolateCrestSpokesQuartic();

    // Superclass method to update the pipeline
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

private:


    int m_InterpolationLevel;
    vtkSRep *m_Input;
    vtkSmartPointer<vtkSRep> m_SRepOutput;

    vtkIdType m_AtomId;
    vtkIdType m_SpokeType;
    double m_Gamma_t;
    double m_Gamma_theta;
    bool m_UseAllSpokes;
    bool m_CyclicSpokes;

    void InterpolateCrest(unsigned numpoints, vector<vtkSRep::SPOKES_TYPE > spokestype,
                          vtkSmartPointer< vtkInterpolateCurve > medialcrestcurveinterpolator,
                          vtkSmartPointer<vtkCellArray>& interpolatedcellarray, vtkSmartPointer<vtkPoints>& interpolatedpoints,
                          vector< vtkSmartPointer< vtkInterpolateCurve > > spokesinterpolator, vector< vtkSmartPointer< vtkInterpolateCurve > > radiusinterpolator,
                          bool usegamma = false);



    vtkSRep::VectorVNLType GetVectorSegment(vtkSRep::VectorVNLType input, int pos);

};

#endif // VTKSREPINTERPOLATEMEDIALSHEET_H
