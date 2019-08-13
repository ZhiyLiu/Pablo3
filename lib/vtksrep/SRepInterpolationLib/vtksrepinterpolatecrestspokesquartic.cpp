#include "vtksrepinterpolatecrestspokesquartic.h"


#include "vtkSmartPointer.h"


#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkQuad.h"

#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vnl/vnl_least_squares_function.h"
#include "vtkAppendPolyData.h"

#include "minimizecurvaturefunction.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

#include "vtkObjectFactory.h"
vtkStandardNewMacro(vtkSRepInterpolateCrestSpokesQuartic);



vtkSRepInterpolateCrestSpokesQuartic::vtkSRepInterpolateCrestSpokesQuartic()
{
    m_InterpolationLevel = 0;
    m_SRepOutput = 0;

    m_AtomId = -1;
    m_SpokeType = -1;
    m_Gamma_t = 1;
    m_Gamma_theta = 1;
    m_UseAllSpokes = false;
    m_CyclicSpokes = false;
}

vtkSRepInterpolateCrestSpokesQuartic::~vtkSRepInterpolateCrestSpokesQuartic()
{
}


// Superclass method to update the pipeline
int vtkSRepInterpolateCrestSpokesQuartic::RequestData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector){



    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkSRep *input = dynamic_cast<vtkSRep*>(vtkSRep::SafeDownCast(inInfo->Get(vtkSRep::DATA_OBJECT())));
    m_Input = input;
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    m_SRepOutput = vtkSmartPointer<vtkSRep>::New();

    vtkSmartPointer<vtkCellArray> interpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> interpolatedpoints = vtkSmartPointer<vtkPoints>::New();


    vtkSRep::VectorIdsType crestids = input->GetCrestMedialAtomsIds();

    if(crestids.size() == 0){
        cout<<"ERROR: No crest atoms found!!!"<<endl;
        return 0;
    }

    vector< vtkSRep::SPOKES_TYPE > spokestype;

    if(m_UseAllSpokes){
        for(int i = 0; i < input->GetNumberOfSpokes(); i++){
            spokestype.push_back((vtkSRep::SPOKES_TYPE)i);
        }
    }else{
        if(m_SpokeType == vtkSRep::BOTTOM_SPOKE){
            spokestype.push_back(vtkSRep::CREST_SPOKE);
            spokestype.push_back(vtkSRep::BOTTOM_SPOKE);

        }else if(m_SpokeType == vtkSRep::TOP_SPOKE){
            spokestype.push_back(vtkSRep::CREST_SPOKE);
            spokestype.push_back(vtkSRep::TOP_SPOKE);

        }else{
            spokestype.push_back(vtkSRep::BOTTOM_SPOKE);
            spokestype.push_back(vtkSRep::CREST_SPOKE);
            spokestype.push_back(vtkSRep::TOP_SPOKE);
        }
    }

    vtkSRep::VectorVNLType crestpositions = input->GetCrestMedialAtoms(crestids);
    vtkSRep::VectorVNLType crestderivatives = input->GetCrestMedialAtomsDerivatives(crestids);
    vtkSRep::VectorVNLType crestnormals = input->GetMedialSheetNormals(crestids);
    vtkSRep::VectorVNLType crestderivativesnormals;

    for(unsigned i = 0; i < crestderivatives.size(); i++){

        VNLType vect = vnl_cross_3d(crestnormals[i], crestderivatives[i]);
        vect = vnl_cross_3d(vect, crestnormals[i]);

        crestderivativesnormals.push_back(vect);
    }



    vector< vtkSmartPointer< vtkInterpolateCurve > > spokesinterpolator;
    vector< vtkSmartPointer< vtkInterpolateCurve > > radiusinterpolator;

    int atomPos = -1;
    unsigned numpoints = crestids.size();
    if(m_AtomId!=-1){
        for(unsigned i = 0; i < crestids.size(); i++){
            if(m_AtomId == crestids[i]){
                atomPos = i;
                numpoints = 2;
            }
        }
        if(atomPos == -1){
            if(this->GetDebug()){
                cout<<"ERROR: Crest interpolation called for an atom that doesn't belong to the crest or end atom!!!"<<endl;
                return 0;
            }else{
                return 1;
            }
        }
    }


    bool cyclic = true;

    if(atomPos != -1){
        cyclic = false;
        crestpositions = GetVectorSegment(crestpositions, atomPos);
        crestderivativesnormals = GetVectorSegment(crestderivativesnormals, atomPos);
    }

    vtkSmartPointer< vtkInterpolateCurve > medialcrestcurveinterpolator = vtkSmartPointer< vtkInterpolateCurve >::New();
    medialcrestcurveinterpolator->SetInterpolationLevel(m_InterpolationLevel);
    medialcrestcurveinterpolator->SetCurvePoints(crestpositions);
    medialcrestcurveinterpolator->SetCurveDerivatives(crestderivativesnormals);
    medialcrestcurveinterpolator->SetCyclic(cyclic);
    medialcrestcurveinterpolator->Update();


    for(unsigned i = 0; i < spokestype.size(); i++){

        vtkSRep::VectorVNLType spokes = input->GetSpokes(crestids, spokestype[i], true);
        vtkSRep::VectorVNLType spokesderivatives = input->GetCrestUDerivatives(crestids, false, spokestype[i]);
        vtkSRep::VectorVNLType spokesradius;

        for(unsigned j = 0; j < crestids.size(); j++){
            vtkSRep::VNLType radius(1);
            radius[0] = input->GetSpokeRadius(crestids[j], spokestype[i]);
            spokesradius.push_back(radius);
        }

        vtkSRep::VectorVNLType spokesradiusderivatives = input->GetCrestMedialAtomsDerivatives(spokesradius);

        if(atomPos != -1){
            spokes = GetVectorSegment(spokes,atomPos);
            spokesderivatives = GetVectorSegment(spokesderivatives,atomPos);
            spokesradius = GetVectorSegment(spokesradius,atomPos);
            spokesradiusderivatives = GetVectorSegment(spokesradiusderivatives,atomPos);
        }else{
            spokes.push_back(spokes[0]);
            spokesderivatives.push_back(spokesderivatives[0]);
            spokesradius.push_back(spokesradius[0]);
            spokesradiusderivatives.push_back(spokesradiusderivatives[0]);
        }

        vtkSmartPointer< vtkInterpolateCurve > spokeinterpolate = vtkSmartPointer< vtkInterpolateCurve >::New();
        spokeinterpolate->SetInterpolationLevel(m_InterpolationLevel);
        spokeinterpolate->SetCurvePoints(spokes);
        spokeinterpolate->SetCurveDerivatives(spokesderivatives);
        spokeinterpolate->Update();

        spokesinterpolator.push_back(spokeinterpolate);

        vtkSmartPointer< vtkInterpolateCurve > radiusinterpolate = vtkSmartPointer< vtkInterpolateCurve >::New();
        radiusinterpolate->SetInterpolationLevel(m_InterpolationLevel);
        radiusinterpolate->SetCurvePoints(spokesradius);
        radiusinterpolate->SetCurveDerivatives(spokesradiusderivatives);
        radiusinterpolate->Update();

        radiusinterpolator.push_back(radiusinterpolate);

    }




    if(atomPos != -1){
        InterpolateCrest(numpoints, spokestype, medialcrestcurveinterpolator, interpolatedcellarray, interpolatedpoints
                         ,spokesinterpolator, radiusinterpolator, true);
    }else{
        InterpolateCrest(numpoints, spokestype, medialcrestcurveinterpolator, interpolatedcellarray, interpolatedpoints
                         ,spokesinterpolator, radiusinterpolator);
    }


    //output->DeepCopy(appendpoly->GetOutput());
    output->SetPoints(interpolatedpoints);
    //output->SetLines(interpolatedcellarray);
    output->SetPolys(interpolatedcellarray);



    return 1;

}

vtkSRep::VectorVNLType vtkSRepInterpolateCrestSpokesQuartic::GetVectorSegment(vtkSRep::VectorVNLType input, int pos){
    vtkSRep::VectorVNLType output;

    if(pos == 0){
        output.push_back(input[input.size()-1]);
        int lim = 2;
        if(lim > (int)input.size()){
            lim = input.size();
        }
        for(int i = 0; i < lim; i++){
            output.push_back(input[i]);
        }
    }else if(pos == ((int)input.size()) - 1){
        for(int i = pos - 1; i < (int)input.size(); i++){
            output.push_back(input[i]);
        }
        output.push_back(input[0]);
    }else{
        int lim = pos + 2;
        if(lim > (int)input.size()){
            lim = input.size();
        }
        for(int i = pos - 1; i < lim; i++){
            output.push_back(input[i]);
        }
    }

    return output;
}

void vtkSRepInterpolateCrestSpokesQuartic::InterpolateCrest(unsigned numpoints, vector<vtkSRep::SPOKES_TYPE > spokestype,
                                                            vtkSmartPointer< vtkInterpolateCurve > medialcrestcurveinterpolator,
                                                            vtkSmartPointer<vtkCellArray>& interpolatedcellarray, vtkSmartPointer<vtkPoints>& interpolatedpoints,
                                                            vector< vtkSmartPointer< vtkInterpolateCurve > > spokesinterpolator, vector< vtkSmartPointer< vtkInterpolateCurve > > radiusinterpolator,
                                                            bool usegamma){

    vtkIdType srepid0 = -1;
    vtkIdType srepidprev = -1;

    vtkSRep::VectorVNLType inplanespokes(spokestype.size());
    vtkSRep::VectorVNLType inplanederivativesspokes(spokestype.size());


    vtkSmartPointer<vtkCellArray> srepoutinterpolatedcellarray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> srepoutinterpolatedpoints = vtkSmartPointer<vtkPoints>::New();
    vtkSRep::RadiusVectorType srepoutallradius;
    vtkSRep::SpokesVectorType srepoutallspokes;

    if(m_Gamma_t > 1){
        m_Gamma_t = 1;
    }
    if(m_Gamma_theta > 1){
        m_Gamma_theta = 1;
    }

    double step = pow((double)2, (double)m_InterpolationLevel);
    double stepsize = 1/step;

    for(unsigned i = 0; i < numpoints; i++){
        double t0 = 0;
        double tend = 1;
        if(usegamma){
            if(i == 0){
                t0 = 1.0 - m_Gamma_t;
                tend = 1;
            }else if(i == numpoints - 1){
                t0 = 0;
                tend = m_Gamma_t;
            }
        }
        vtkSRep::VectorSRepIdsType pointsIds;
        int pointidindex = -1;

        for(double t = t0; t <= tend; t+=stepsize){

            VNLType p0 = medialcrestcurveinterpolator->EvaluateFunction(i, t);
            //vtkIdType id0 = interpolatedpoints->InsertNextPoint(p0[0], p0[1], p0[2]);

            if(srepid0 != -1){
                srepidprev = srepid0;
            }
            srepid0 = srepoutinterpolatedpoints->InsertNextPoint(p0[0], p0[1], p0[2]);

            if(srepidprev != -1){
                vtkSmartPointer<vtkLine> srepline = vtkSmartPointer<vtkLine>::New();
                srepline->GetPointIds()->SetId(0, srepidprev);
                srepline->GetPointIds()->SetId(1, srepid0);
                srepoutinterpolatedcellarray->InsertNextCell(srepline);
            }

            for(unsigned j = 0; j < spokestype.size(); j++){

                VNLType spoke = spokesinterpolator[j]->EvaluateFunction(i, t);
                VNLType radius = radiusinterpolator[j]->EvaluateFunction(i, t);

                VNLType p1 = p0 + spoke*radius[0];

                inplanespokes[j] = p1;

                /*vtkIdType id1 = interpolatedpoints->InsertNextPoint(p1[0], p1[1], p1[2]);

                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                line->GetPointIds()->SetId(0, id0);
                line->GetPointIds()->SetId(1, id1);
                interpolatedcellarray->InsertNextCell(line);*/
            }


            //Calculate the derivatives in "plane" of the three spokes
            for(unsigned j = 0; j < spokestype.size(); j++){
                if( j == 0){
                    VNLType t = inplanespokes[j];
                    if(m_CyclicSpokes){
                        t = (inplanespokes[j+1] - inplanespokes[spokestype.size()-1])/2.0;
                    }else{
                        t.fill(0);
                    }
                    inplanederivativesspokes[j] = t;
                }else if(j == spokestype.size() - 1){                    
                    VNLType t = inplanespokes[j];
                    if(m_CyclicSpokes){
                        t = (inplanespokes[0] - inplanespokes[j-1])/2.0;
                    }else{
                        t.fill(0);
                    }
                    inplanederivativesspokes[j] = t;
                }else{
                    inplanederivativesspokes[j] = (inplanespokes[j+1] - inplanespokes[j-1])/2.0;
                }
            }

            //fit the crest curve for the plane
            vtkSmartPointer< vtkInterpolateCurve > inplaneinterpolate = vtkSmartPointer< vtkInterpolateCurve >::New();
            inplaneinterpolate->SetInterpolationLevel(m_InterpolationLevel);
            inplaneinterpolate->SetCurvePoints(inplanespokes);
            inplaneinterpolate->SetCurveDerivatives(inplanederivativesspokes);
            inplaneinterpolate->SetCyclic(m_CyclicSpokes);
            inplaneinterpolate->Update();


            vtkSRep::VectorVNLType srepoutspokes;
            vtkSRep::VectorDoubleType srepoutradius;


            pointsIds.push_back(vtkSRep::VectorIdsType());
            pointidindex++;

            int limit = spokestype.size()-1;
            if(m_CyclicSpokes){
                limit = spokestype.size();
            }

            for(int j = 0; j < limit; j++){
                double theta0 = 0;
                double thetaend = 1;
                if(j == 0){
                    theta0 = 1.0 - m_Gamma_theta;
                    thetaend = 1;
                }else if(j == ((int)spokestype.size()) - 2){
                    theta0 = 0;
                    thetaend = m_Gamma_theta;
                }                

                for(double theta = theta0; theta <= thetaend; theta+=stepsize){

                    VNLType p1 = inplaneinterpolate->EvaluateFunction(j, theta);

                    //vtkIdType id1 = interpolatedpoints->InsertNextPoint(p1[0], p1[1], p1[2]);
                    //vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    //line->GetPointIds()->SetId(0, id0);
                    //line->GetPointIds()->SetId(1, id1);
                    //interpolatedcellarray->InsertNextCell(line);

                    if(theta != 0 || j != 0){
                        vtkSRep::VNLType srepoutspoke = p1 - p0;
                        srepoutradius.push_back(srepoutspoke.magnitude());
                        srepoutspoke = srepoutspoke.normalize();

                        srepoutspokes.push_back(srepoutspoke);
                    }



                    pointsIds[pointidindex].push_back(interpolatedpoints->InsertNextPoint(p1[0], p1[1], p1[2]));
                }                
            }

            srepoutallspokes.push_back(srepoutspokes);
            srepoutallradius.push_back(srepoutradius);

            srepoutspokes.clear();
            srepoutradius.clear();

        }

        for(unsigned j = 0; j < pointsIds.size() - 1; j++){
            for(unsigned k = 0; k < pointsIds[j].size() - 1; k++){
                vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                quad->GetPointIds()->SetId(0, pointsIds[j][k]);
                quad->GetPointIds()->SetId(1, pointsIds[j+1][k]);
                quad->GetPointIds()->SetId(2, pointsIds[j+1][k+1]);
                quad->GetPointIds()->SetId(3, pointsIds[j][k+1]);

                //quad->Print(cout);

                interpolatedcellarray->InsertNextCell(quad);
            }
        }
    }

    m_SRepOutput->SetPoints(srepoutinterpolatedpoints);
    m_SRepOutput->SetLines(srepoutinterpolatedcellarray);
    m_SRepOutput->SetAllSpokesRadius(srepoutallradius);
    m_SRepOutput->SetAllSpokes(srepoutallspokes);
}
