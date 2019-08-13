#ifndef optimizeCortexSpokes_H
#define optimizeCortexSpokes_H

#include "vnl/vnl_least_squares_function.h"

#include <vtkDataArray.h>
#include <vtksrep.h>
#include <vtkPointLocator.h>
#include <vtkSmartPointer.h>
#include "vtksrepinterpolatemedialspokeshermite.h"

class optimizeCortexSpokes : public vnl_least_squares_function
{
public:
    optimizeCortexSpokes();
    ~optimizeCortexSpokes();

    virtual void f(vnl_vector< double > const &x, vnl_vector< double > &fx);

    void SetPointLocator0(vtkSmartPointer<vtkPointLocator> locator){
        m_Locator0 = locator;
    }    

    void SetPointLocator1(vtkSmartPointer<vtkPointLocator> locator){
        m_Locator1 = locator;
    }

    void SetAtomId(vtkIdType atomId){
        m_AtomId = atomId;
    }

    void SetSrep(vtkSRep* srep){
        m_Srep = srep;
    }

private:


    vtkSmartPointer<vtkPointLocator> m_Locator0;
    vtkSmartPointer<vtkPointLocator> m_Locator1;

    vtkIdType m_AtomId;
    vtkSRep* m_Srep;

    vtkSmartPointer< vtkSRepInterpolateMedialSpokesHermite > m_Medialspokesinterpolator;

};

#endif // optimizeCortexSpokes_H
