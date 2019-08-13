#include "vtkoptimizeSRepSampling.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkQuad.h"



#include "vtkInformationVector.h"
#include "vtkInformation.h"





//#include "vnl/algo/vnl_conjugate_gradient.h"
#include "vnl/algo/vnl_levenberg_marquardt.h"
#include "vnl/vnl_least_squares_cost_function.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkOptimizeSRepSampling);

#include <map>

#include "optimizeSRepSampling.h"

vtkOptimizeSRepSampling::vtkOptimizeSRepSampling()
{


}

vtkOptimizeSRepSampling::~vtkOptimizeSRepSampling()
{
}


// Superclass method to update the pipeline
void vtkOptimizeSRepSampling::Update(){

    m_SRep = vtkSmartPointer<vtkSRep>::New();
    vtkSmartPointer< vtkPoints > hubpos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSRep::VectorSRepIdsType pointsIds;

    vtkSRep::RadiusVectorType allradius;
    vtkSRep::SpokesVectorType allspokes;

    vector<double> steps;
    double numstep = (m_Max - m_Min)/((double)m_NumP);

    double step = 0.0;
    for(step = m_Min; step < m_Max; step+=numstep){//40*2
        steps.push_back(step);
    }    
    steps.insert(steps.end(), steps.rbegin(), steps.rend());
    steps.insert(steps.begin() + m_NumP, 1.0);

    /*steps.clear();
    numstep = 1.0/m_NumP;
    logsum = 0;
    for(double i = 0; i < 1; i+=numstep){
        steps.push_back( i);
        logsum += i;
    }
    steps.insert(steps.end(), steps.rbegin(), steps.rend());*/

    for(int xn = 0; xn < steps.size(); xn++){
        pointsIds.push_back(vtkSRep::VectorIdsType());
        for(int yn = 0; yn < steps.size(); yn++){

            double x = steps[xn];
            double y = steps[yn];

            double X = log(x)/log(m_Min);
            if(xn >= m_NumP){
                X = -X;
            }
            double Y = log(y)/log(m_Min);
            if(yn >= m_NumP){
                Y = -Y;
            }
            double temp = X * sqrt(1 - pow(Y, 2.0)/2.0);

            Y = Y * sqrt( 1 - pow(X,2)/2)*4;
            X = temp*4;
            //X *= 4;
            //Y *= 4;

            double point[3];
            //stereographic
            point[0] = (2*X/(1+pow(X,2)+pow(Y,2)));
            point[1] = (2*Y/(1+pow(X,2)+pow(Y,2)));
            point[2] = ((-1+pow(X,2)+pow(Y,2))/(1+pow(X,2)+pow(Y,2)));
            /*point[0] = X;
            point[1] = Y;
            point[2] = 0;*/


            vnl_vector<double> midpoint(point, 3);
            midpoint*=100;

            pointsIds[(int)xn].push_back(hubpos->InsertNextPoint(midpoint[0], midpoint[1], midpoint[2]));

            vtkSRep::VNLType s(3);
            s.fill(0.1);
            vtkSRep::VNLType s1(3);
            s1.fill(-0.1);

            vtkSRep::VectorDoubleType radius;
            radius.push_back(s.magnitude());
            radius.push_back(s1.magnitude());

            vtkSRep::VectorVNLType vnlspokes;

            vnlspokes.push_back(s.normalize());
            vnlspokes.push_back(s1.normalize());

            bool insertendspoke = false;
            if(xn == 0 || xn == steps.size() || yn == 0 || yn == steps.size()){
                insertendspoke = true;
            }

            if(insertendspoke){
                vtkSRep::VNLType send(3);
                send.fill(0.1);

                radius.push_back(send.magnitude());
                vnlspokes.push_back(send.normalize());
            }

            allspokes.push_back(vnlspokes);
            allradius.push_back(radius);

            /*if(yn < steps.size()){
                y += steps[yn];
            }*/
        }
        /*if(xn < steps.size()){
            x += steps[xn];
        }*/
    }


    for(unsigned i = 0; i < pointsIds.size() - 1; i++){
         for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

             vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
             quad->GetPointIds()->SetId(0, pointsIds[i][j]);
             quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
             quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
             quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

             //quad->Print(cout);
             cellarray->InsertNextCell(quad);

         }
     }

    m_SRep->SetPoints(hubpos);
    m_SRep->SetPolys(cellarray);
    //srepfig->BuildLinks();
    m_SRep->SetAllSpokes(allspokes);
    m_SRep->SetAllSpokesRadius(allradius);
    //srepfig->GetPointData()->SetNormals(normalarray);


    if(m_SRep->GetNumberOfPoints() == (int)m_RegisteredPoints.size()){

        for(unsigned i = 0; i < m_RegisteredPoints.size(); i++){
            VNLType p = m_RegisteredPoints[i];

            //p[0] = i/m_SRep->GetNumColumns() - (m_SRep->GetNumColumns()-1)/2.0;
            //p[1] = i%m_SRep->GetNumColumns() - (m_SRep->GetNumRows()-1)/2.0;

            double xstep = p[0];
            double ystep = p[1];


            xstep = 2*fabs(p[0])/(m_SRep->GetNumColumns()-1) * (m_Min - m_Max) + m_Max;
            ystep = 2*fabs(p[1])/(m_SRep->GetNumRows()-1) * (m_Min - m_Max) + m_Max;
            //m_SRep->GetPoints()->SetPoint(i, xstep, ystep, 0);

            if(xstep < m_Min){
                xstep = m_Min;
            }

            double X = log(xstep)/log(m_Min);

            if(p[0] > 0){
                X = -X;
            }


            if(ystep < m_Min){
                ystep = m_Min;
            }
            double Y = log(ystep)/log(m_Min);
            if(p[1] > 0){
                Y = -Y;
            }



            double temp = X * sqrt(1 - pow(Y, 2.0)/2.0);

            if(isnan(temp)){
                cout<<"isnan"<<endl;
            }

            Y = Y * sqrt( 1 - pow(X,2)/2)*4;
            X = temp*4;

            double point[3];
            //stereographic
            point[0] = (2*X/(1+pow(X,2)+pow(Y,2)));
            point[1] = (2*Y/(1+pow(X,2)+pow(Y,2)));
            point[2] = ((-1+pow(X,2)+pow(Y,2))/(1+pow(X,2)+pow(Y,2)));



            //cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;

            /*point[0] = p[0];
            point[1] = p[1];
            point[2] = 0;*/

            //if(xn <= m_NumP){
            //m_SRep->GetPoints()->SetPoint(curri, xstep*10, ystep*10, 0);
            m_SRep->GetPoints()->SetPoint(i, point[0], point[1], point[2]);
        }

        /*for(unsigned xn = 0; xn < m_SRep->GetNumColumns(); xn++){
            pointsIds.push_back(vtkSRep::VectorIdsType());
            for(unsigned yn = 0; yn < m_SRep->GetNumColumns(); yn++){

                int curri = xn*m_SRep->GetNumColumns() + yn;
                double p[3];
                m_SRep->GetPoint(curri, p);

                double x = p[0];
                double y = p[1];

                double X = log(x)/log(m_Min);
                if(xn >= m_NumP){
                    X = -X;
                }
                double Y = log(y)/log(m_Min);
                if(yn >= m_NumP){
                    Y = -Y;
                }


                double temp = X * sqrt(1 - pow(Y, 2.0)/2.0);

                Y = Y * sqrt( 1 - pow(X,2)/2)*4;
                X = temp*4;

                double point[3];
                //stereographic
                point[0] = (2*X/(1+pow(X,2)+pow(Y,2)));
                point[1] = (2*Y/(1+pow(X,2)+pow(Y,2)));
                point[2] = ((-1+pow(X,2)+pow(Y,2))/(1+pow(X,2)+pow(Y,2)));

                //if(xn <= m_NumP){
                //m_SRep->GetPoints()->SetPoint(curri, xstep*10, ystep*10, 0);
                m_SRep->GetPoints()->SetPoint(curri, point[0], point[1], point[2]);


            }
        }*/
    }
}

/*
    vector<double> steps;
    double logsum = 0;

    double numstep = (m_Max - m_Min)/((double)m_NumP);

    for(double step = m_Min; step <= m_Max; step+=numstep){//40*2
        //cout<<step<<endl;
        double temp = fabs(log10(step)/log10(m_Min));
        steps.push_back(temp);
        logsum += temp;
    }
    steps.insert(steps.end(), steps.rbegin(), steps.rend());

    //steps.clear();
    //numstep = 1.0/m_NumP;
    //logsum = 0;
    //for(double i = 0; i < 1; i+=numstep){
    //    steps.push_back( i);
    //    logsum += i;
    //}
    //steps.insert(steps.end(), steps.rbegin(), steps.rend());

    double x = -logsum;
    for(unsigned xn = 0; xn <= steps.size(); xn++){
        double y = -logsum;
        pointsIds.push_back(vtkSRep::VectorIdsType());
        for(unsigned yn = 0; yn <= steps.size(); yn++){

            double X = x/logsum;
            double Y = y/logsum;
            double temp = X * sqrt(1 - pow(Y, 2.0)/2.0);

            Y = Y * sqrt( 1 - pow(X,2)/2)*4;
            X = temp*4;
            //X *= 4;
            //Y *= 4;

            double point[3];
            //stereographic
            point[0] = (2*X/(1+pow(X,2)+pow(Y,2)));
            point[1] = (2*Y/(1+pow(X,2)+pow(Y,2)));
            point[2] = ((-1+pow(X,2)+pow(Y,2))/(1+pow(X,2)+pow(Y,2)));
            //point[0] = X;
            //point[1] = Y;
            //point[2] = 0;


            vnl_vector<double> midpoint(point, 3);
            midpoint*=100;

            pointsIds[(int)xn].push_back(hubpos->InsertNextPoint(midpoint[0], midpoint[1], midpoint[2]));

            vtkSRep::VNLType s(3);
            s.fill(0.1);
            vtkSRep::VNLType s1(3);
            s1.fill(-0.1);

            vtkSRep::VectorDoubleType radius;
            radius.push_back(s.magnitude());
            radius.push_back(s1.magnitude());

            vtkSRep::VectorVNLType vnlspokes;

            vnlspokes.push_back(s.normalize());
            vnlspokes.push_back(s1.normalize());

            bool insertendspoke = false;
            if(xn == 0 || xn == steps.size() || yn == 0 || yn == steps.size()){
                insertendspoke = true;
            }

            if(insertendspoke){
                vtkSRep::VNLType send(3);
                send.fill(0.1);

                radius.push_back(send.magnitude());
                vnlspokes.push_back(send.normalize());
            }

            allspokes.push_back(vnlspokes);
            allradius.push_back(radius);

            if(yn < steps.size()){
                y += steps[yn];
            }
        }
        if(xn < steps.size()){
            x += steps[xn];
        }
    }


    for(unsigned i = 0; i < pointsIds.size() - 1; i++){
         for(unsigned j = 0; j < pointsIds[i].size() - 1; j++){

             vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
             quad->GetPointIds()->SetId(0, pointsIds[i][j]);
             quad->GetPointIds()->SetId(1, pointsIds[i+1][j]);
             quad->GetPointIds()->SetId(2, pointsIds[i+1][j+1]);
             quad->GetPointIds()->SetId(3, pointsIds[i][j+1]);

             //quad->Print(cout);
             cellarray->InsertNextCell(quad);

         }
     }

    m_SRep->SetPoints(hubpos);
    m_SRep->SetPolys(cellarray);
    //srepfig->BuildLinks();
    m_SRep->SetAllSpokes(allspokes);
    m_SRep->SetAllSpokesRadius(allradius);
    //srepfig->GetPointData()->SetNormals(normalarray);
*/
