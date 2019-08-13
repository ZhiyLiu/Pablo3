#ifndef OptimizeSRepSampling_H
#define OptimizeSRepSampling_H

//#include "vnl/vnl_cost_function.h"
#include "vnl/vnl_least_squares_function.h"

#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkPointLocator.h>
#include <vtkSmartPointer.h>
#include <vtksrep.h>

class OptimizeSRepSampling : public vnl_least_squares_function
{
public:
    OptimizeSRepSampling(int numcells, int numpoints);
    ~OptimizeSRepSampling();

    virtual void f(vnl_vector< double > const &x, vnl_vector< double > &fx);
    //virtual double f(vnl_vector< double > const &x);
    //virtual void gradf(vnl_vector< double > const &x, vnl_vector< double > &gradient);

    void SetSRep(vtkSRep* srep){
        m_SRep = srep;
    }

    void SetNumP(int nump){
        m_NumP = nump;
    }

private:

    vtkSRep* m_SRep;
    int m_NumP;

};

#endif // OptimizeSRepSampling_H
