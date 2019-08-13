#ifndef VTKSREP_H
#define VTKSREP_H

#include "vtkPolyData.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_inverse.h"
#include <vector>

using namespace std;

#define CURVATURE_ARRAY "CURVATURE_ARRAY"

class vtkSRep : public vtkPolyData
{
public:
    static vtkSRep *New();

    typedef vector< double > VectorDoubleType;
    typedef vector< VectorDoubleType > RadiusVectorType;

    typedef vnl_vector<double> VNLType;
    typedef vector< VNLType > VectorVNLType;
    typedef vector< VectorVNLType > SpokesVectorType;

    typedef vector< vtkIdType >     VectorIdsType;
    typedef vector< VectorIdsType > VectorSRepIdsType;

    enum SPOKES_TYPE { TOP_SPOKE = 0, BOTTOM_SPOKE = 1, CREST_SPOKE = 2 };    

    typedef vnl_matrix<double> VNLMatrixType;

    void SetAllSpokesRadius(RadiusVectorType radius){
        m_SpokesRadius = radius;
    }

    RadiusVectorType GetAllRadius(){
        return m_SpokesRadius;
    }

    VectorDoubleType GetSpokesRadius(vtkIdType hubid){
        if(hubid > (int)m_SpokesRadius.size()){
            cout<<"hubid does not exists"<<endl;

        }else{
            return m_SpokesRadius[hubid];
        }
        return VectorDoubleType();
    }

    double GetSpokeRadius(vtkIdType hubid, vtkIdType spokepos){
        if(hubid > (int)m_SpokesRadius.size()){
            cout<<"hubid does not exists"<<endl;

        }else{
            return m_SpokesRadius[hubid][spokepos];
        }
        return -1;
    }

    void SetAllSpokes(SpokesVectorType spokes){
        m_Spokes = spokes;
    }

    SpokesVectorType GetAllSpokes(){
        return m_Spokes;
    }

    VNLType GetTopSpoke(vtkIdType hubid){
        if(hubid > (int)m_Spokes.size()){
            cout<<"hubid does not exists"<<endl;
            return VNLType(3);
        }else{
            return m_Spokes[hubid][TOP_SPOKE] * m_SpokesRadius[hubid][TOP_SPOKE];
        }

    }

    VNLType GetBottomSpoke(vtkIdType hubid){
        if(hubid > (int)m_Spokes.size()){
            cout<<"hubid does not exists"<<endl;
            return VNLType(3);
        }else{
            return m_Spokes[hubid][BOTTOM_SPOKE] * m_SpokesRadius[hubid][BOTTOM_SPOKE];
        }

    }

    VectorVNLType GetSpokes(vtkIdType hubid){
        if(hubid > (int)m_Spokes.size()){
            cout<<"hubid does not exists"<<endl;
        }else{
            return m_Spokes[hubid];
        }
        return VectorVNLType();
    }
    /* \brief Given a vector of hub ids and the spoke position
              returns a list containing the hub spoke at the spokeposition
              (each hub contains 1..N spokes)
    */
    VectorVNLType GetSpokes(VectorIdsType vecthubids, vtkIdType spokepos, bool unit = false);

    VNLType GetSpoke(vtkIdType hubid, vtkIdType spokepos, bool unit = false){
        if(hubid > (int)m_Spokes.size()){
            cout<<"hubid does not exists"<<endl;
        }else{
            VNLType spoke = m_Spokes[hubid][spokepos]*m_SpokesRadius[hubid][spokepos];
            if(unit){
                return spoke.normalize();
            }
            return spoke;
        }
        return VNLType(0);
    }

    void SetSpoke(vtkIdType hubid, vtkIdType spokepos, VNLType spoke){
        if(hubid > (int)m_Spokes.size()){
            cout<<"hubid does not exists"<<endl;
        }else{
            m_SpokesRadius[hubid][spokepos] = spoke.magnitude();
            m_Spokes[hubid][spokepos] = spoke.normalize();
        }
    }
    void SetSpokes(vtkIdType hubid, VectorVNLType spokes){
        m_Spokes[hubid] = spokes;
    }

    // \brief given an id of a hub returns the e3 vector i.e. cross(topSpoke, bottomSpoke)
    VNLType GetMedialSheetVecte3(vtkIdType hubid);
    // \brief given an id of a hub returns the normal to the sheet
    VNLType GetMedialSheetNormal(vtkIdType hubid);
    // \brief given an id of a hub returns the normal to the sheet
    VectorVNLType GetMedialSheetNormals(VectorIdsType hubids);

    VectorVNLType GetDerivatives(vtkIdType pointid);
    //VectorVNLType GetDerivatives(vtkIdType rowid, vtkIdType colid);
    VectorVNLType GetDerivativesSpoke(vtkIdType pointid, vtkIdType spokepos, bool unit = false);
    VectorVNLType GetDerivativesSpokeRadius(vtkIdType pointid, vtkIdType spokepos);

    /*void SetGridTopolgyIds(VectorSRepIdsType gridids){
        m_GridTopolgyIds = gridids;
    }*/
    /*VectorSRepIdsType GetGridTopolgyIds(){
        return m_GridTopolgyIds;
    }*/

    // \brief returns a list of vnl_vectors that contain the position of every crest medial atom
    VectorVNLType GetCrestMedialAtoms();
    // \brief returns a list of vnl_vectors that contain the position of every crest medial atom
    // \param crestids the list of the ids of crestmedialAtoms
    VectorVNLType GetCrestMedialAtoms(VectorIdsType crestids);
    // \brief returns a list of vnl_vectors that contain the ids of every crest medial atom
    VectorIdsType GetCrestMedialAtomsIds();
    // \brief returns a list of vnl_vectors that contain the ids of every crest medial atom on the given slab side
    //        side [0-3]
    VectorIdsType GetCrestMedialAtomsIds(int side);
    // \brief returns a list of vnl_vectors that contain the medial atoms finite derivatives (tangent vectors to the points)
    // \params bool unit = false
    VectorVNLType GetCrestMedialAtomsDerivatives(bool unit = false);
    VectorVNLType GetCrestMedialAtomsDerivatives(VectorIdsType crestids, bool unit = false);
    VectorVNLType GetCrestMedialAtomsDerivatives(VectorVNLType crestpositions, bool unit = false);

    // \brief returns a list that contains the finite derivatives of the crest spoke
    // \param bool unit as default the list return the value of the derivative without making it unit
    VectorVNLType GetCrestUDerivatives(bool unit = false, vtkIdType spokepos = CREST_SPOKE);
    VectorVNLType GetCrestUDerivatives(VectorIdsType crestids, bool unit = false, vtkIdType spokepos = CREST_SPOKE);
    VectorVNLType GetCrestUDerivatives(VectorVNLType crestvectors, bool unit = false);

    VNLMatrixType GetSRadMatrix(vtkIdType hubid, vtkIdType spokepos);

    VectorIdsType GetInternalMedialAtomIds(vtkIdType startid, vtkIdType nextincellid);


    //global number of spokes
    //Set it when using quasitubes
    void SetNumberOfSpokes(int numspokes){
        m_NumberOfSpokes = numspokes;
    }
    int GetNumberOfSpokes(){
        return m_NumberOfSpokes;
    }


    double GetCellArea(vtkIdType cellid);
    double GetCellAreaAverage();

    int GetNumRows();
    int GetNumColumns();
protected:
    vtkSRep();

    ~vtkSRep();

private:

    SpokesVectorType m_Spokes;
    RadiusVectorType m_SpokesRadius;

    //VectorSRepIdsType m_GridTopolgyIds;
    int m_NumberOfSpokes;

    /* \brief Calculates the next atom id that corresponds to a crest medial atom
       \params vtkIdType id  the current atom id that is to evaluate
               this pointid is use to recover all the cells that this point belong to
       \params vtkIdType inquadid corresponds to the current inquadid for the given id
               this is use to locate the correct cell
               if no cell is found to have the inquadid for the given id then
               inquadid is increased. This happens at the corners of the slab
        \params vectids returns all the ids of the points found
    */
    void GetCrestMedialAtomsIds(vtkIdType pointid, vtkIdType inquadid, VectorIdsType& vectids );

    void GetInternalMedialAtomIds(vtkIdType startid, vtkIdType nextincellid, VectorIdsType& vectids );

    int m_NumRows;
    int m_NumColumns;

};

#endif // VTKSREP_H
