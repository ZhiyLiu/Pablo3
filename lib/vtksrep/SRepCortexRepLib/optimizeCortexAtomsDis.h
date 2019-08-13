#ifndef optimizeCortexAtomsDis_H
#define optimizeCortexAtomsDis_H

#include "vnl/vnl_least_squares_function.h"

#include <vtkDataArray.h>
#include <vtksrep.h>
#include <vtkPointLocator.h>
#include <vtkSmartPointer.h>

class optimizeCortexAtomsDis : public vnl_least_squares_function
{
public:
    optimizeCortexAtomsDis();
    ~optimizeCortexAtomsDis();

    virtual void f(vnl_vector< double > const &x, vnl_vector< double > &fx);

    void SetPointLocator(vtkSmartPointer<vtkPointLocator> locator){
        m_Locator = locator;
    }    

    void SetSRep(vtkSRep* srep){
        m_SRep = srep;
    }

    void SetId(vtkIdType id){
        m_Id = id;
    }

private:


    vtkSmartPointer<vtkPointLocator> m_Locator;
    int m_NPoints;
    vtkSRep* m_SRep;
    vtkIdType m_Id;

};

#endif // optimizeCortexAtomsDis_H
