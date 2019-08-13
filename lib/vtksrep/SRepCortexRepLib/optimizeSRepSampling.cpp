#include "optimizeSRepSampling.h"

#include "iostream"

#include "vtkIdList.h"
#include "vtkCell.h"

using namespace std;

OptimizeSRepSampling::OptimizeSRepSampling(int numcells, int numpoints)
    : vnl_least_squares_function(2, numcells, no_gradient)
{
    m_SRep = 0;
    m_NumP = numpoints;
}


OptimizeSRepSampling::~OptimizeSRepSampling()
{
}

void OptimizeSRepSampling::f(vnl_vector< double > const &x, vnl_vector< double > &fx){
    if(m_SRep){
        double min = x[0];
        double max = x[1];

        cout<<x<<endl;

        if(min > 0 && max > 0){

            double numstep = (max - min)/((double)m_NumP);
            vector<double> steps;
            double logsum = 0;



            for(double step = min; step <= max; step+=numstep){
                double temp = fabs(log10(step)/log10(min));
                steps.push_back(temp);
                logsum += temp;
            }
            steps.insert(steps.end(), steps.rbegin(), steps.rend());

            double x = -logsum;
            for(unsigned xn = 0; xn <= steps.size(); xn++){
                double y = -logsum;

                for(unsigned yn = 0; yn <= steps.size(); yn++){

                    double X = x/logsum;
                    double Y = y/logsum;
                    double temp = X * sqrt(1 - pow(Y, 2.0)/2.0);

                    Y = Y * sqrt( 1 - pow(X,2)/2)*4;
                    X = temp*4;

                    double point[3];
                    //stereographic
                    point[0] = (2*X/(1+pow(X,2)+pow(Y,2)));
                    point[1] = (2*Y/(1+pow(X,2)+pow(Y,2)));
                    point[2] = ((-1+pow(X,2)+pow(Y,2))/(1+pow(X,2)+pow(Y,2)));

                    vnl_vector<double> midpoint(point, 3);
                    midpoint*=100;

                    vtkIdType pid = xn*(steps.size()+1) + yn;

                    m_SRep->GetPoints()->SetPoint(pid, midpoint[0], midpoint[1], midpoint[2]);


                    if(yn < steps.size()){
                        y += steps[yn];
                    }
                }
                if(xn < steps.size()){
                    x += steps[xn];
                }
            }
            //double avg = m_SRep->GetCellAreaAverage();
            //double retval = 0;
            for(unsigned i = 0; i < m_SRep->GetNumberOfCells(); i++){
                //double val = m_SRep->GetCellArea(i) - avg;
                fx[i] = m_SRep->GetCellArea(i);
            }

        }else{
            fx.fill(99999);
        }
    }
}
