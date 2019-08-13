#include "itkThinPlateSplineExtended.h"

#include <iostream>
#include <string>

#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "itkPointSet.h"
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkDoubleArray.h>
#include <vtkCenterOfMass.h>
#include <vtkObjectFactory.h>
#include <vtkCurvatures.h>
#include <vtkVector.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPolyDataNormals.h>

#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkMassProperties.h>
#include <vtkCellArray.h>

int main( int argc, char* argv[] ){

    using namespace std;
	typedef double CoordinateRepType;
	typedef itkThinPlateSplineExtended TransformType;
	typedef itk::Point< CoordinateRepType, 3 > PointType;
	typedef std::vector< PointType > PointArrayType;
	typedef TransformType::PointSetType PointSetType;
	typedef PointSetType::Pointer PointSetPointer;
	typedef PointSetType::PointIdentifier PointIdType;

	// PointSetType::Pointer sourceLandMarks = PointSetType::New();
	// PointSetType::Pointer targetLandMarks = PointSetType::New();
	// PointType p1; PointType p2; // same as double p1[3];
	// PointSetType::PointsContainer::Pointer sourceLandMarkContainer
	// 		= sourceLandMarks->GetPoints();
	// PointSetType::PointsContainer::Pointer targetLandMarkContainer
	// 		= targetLandMarks->GetPoints();

	// PointIdType id_s = itk::NumericTraits< PointIdType >::Zero;
	// PointIdType id_t = itk::NumericTraits< PointIdType >::Zero;

	// // Read in the source points set
	// vtkSmartPointer<vtkPolyData> polyData_source; // same as vtkPolyData* polyData_source;
	// vtkSmartPointer<vtkPolyData> polyData_target;

	// vtkSmartPointer<vtkPolyDataReader> reader_source
	// 		= vtkSmartPointer<vtkPolyDataReader>::New();
	// vtkSmartPointer<vtkPolyDataReader> reader_target
	// 		= vtkSmartPointer<vtkPolyDataReader>::New();

	// std::string sourcefilename("/tmp/Slicer-zhiy/forward/0501.vtk");
	// reader_source->SetFileName(sourcefilename.c_str());
	// reader_source->Update();
	// polyData_source = reader_source->GetOutput();
	// cout<<"----sourcefilename: "<<sourcefilename<<endl;

	// std::string targetfilename("/tmp/Slicer-zhiy/forward/0500.vtk");
	// reader_target->SetFileName(targetfilename.c_str());
	// reader_target->Update();
	// polyData_target = reader_target->GetOutput();
	// cout<<"----targetfilename: "<<targetfilename<<endl;

	// // Read in the source points set
	// for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); ++i){
	// 	double p[3];
	// 	polyData_source->GetPoint(i,p);
	// 	p1[0] = p[0];
	// 	p1[1] = p[1];
	// 	p1[2] = p[2];
	// 	sourceLandMarkContainer->InsertElement(id_s, p1);
	// 	id_s++;
	// }

	// // Read in the target points set
	// for(unsigned int i = 0; i < polyData_target->GetNumberOfPoints(); ++i){
	// 	double p[3];
	// 	polyData_target->GetPoint(i,p);
	// 	p2[0] = p[0];
	// 	p2[1] = p[1];
	// 	p2[2] = p[2];
	// 	targetLandMarkContainer->InsertElement(id_t, p2);
	// 	id_t++;
	// }


    // TransformType::Pointer tps = TransformType::New();
	// tps->SetSourceLandmarks(sourceLandMarks);
	// tps->SetTargetLandmarks(targetLandMarks);

	// cout<<"Computing W Matrix... "<<endl;
	// tps->ComputeWMatrix();
    ifstream inFile;
    inFile.open("/tmp/Slicer-zhiy/backward/0499.txt");
    if(!inFile)
    {
        std::cerr << "unable to open file" << std::endl;
        return -1;
    }
    	itkThinPlateSplineExtended::DMatrixType D;
	itkThinPlateSplineExtended::AMatrixType A;
	itkThinPlateSplineExtended::BMatrixType B;
	// first read in the size
	string buffer;
	std::getline(inFile,buffer);
	std::istringstream ss1(buffer);
	int nRows = 0;
	int nCols = 0;
	char tmp;
	ss1 >> nRows >> tmp >> nCols;
	D.set_size(nRows,nCols);
	for(int i = 0; i < nRows; i++) {
		for(int j = 0; j < nCols; j++) {
			string buffer2;
			std::getline(inFile, buffer2, ',');
			double entry = atof(buffer2.c_str());
			D(i,j) = entry;
		}
	}
	buffer.clear();
	std::getline(inFile, buffer);
	std::getline(inFile, buffer);
	for(int i = 0; i < A.rows(); i++) {
		for(int j = 0; j < A.cols(); j++) {
			string buffer2;
			std::getline(inFile, buffer2,',');
			double entry = atof(buffer2.c_str());
			A(i,j) = entry;
		}
	}
	buffer.clear();
	std::getline(inFile, buffer);
	std::getline(inFile, buffer);
	for(int i = 0; i < B.size(); i++) {
		string buffer2;
		std::getline(inFile, buffer2, ',');
		double entry = atof(buffer2.c_str());
		B(i) = entry;
	}

    TransformType::Pointer tps = TransformType::New();
	tps->setDMatrix(D);
	tps->setAMatrix(A);
	tps->setBVector(B);

	PointSetType::Pointer sourceLandMarks = PointSetType::New();
	PointSetType::PointsContainer::Pointer sourceLandMarkContainer = sourceLandMarks->GetPoints();
	vtkSmartPointer<vtkPolyDataReader> reader_source = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkSmartPointer<vtkPolyData> polyData_source = vtkSmartPointer<vtkPolyData>::New();
	reader_source->SetFileName("/tmp/Slicer-zhiy/forward/0500.vtk");
	reader_source->Update();
	polyData_source = reader_source->GetOutput();
	PointIdType id_s = itk::NumericTraits< PointIdType >::Zero;
	PointType p1;
	// Read in the source points set
	for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i += 10){
		double p[3];
		polyData_source->GetPoint(i,p);
		// This is identical to:
		// polydata->GetPoints()->GetPoint(i,p);
		p1[0] = p[0];
		p1[1] = p[1];
		p1[2] = p[2];
		sourceLandMarkContainer->InsertElement(id_s, p1);
		id_s++;
	}
	tps->SetSourceLandmarks(sourceLandMarks);

    vtkSmartPointer<vtkPolyDataReader> srepReader = vtkSmartPointer<vtkPolyDataReader>::New();
    srepReader->SetFileName("/tmp/Slicer-zhiy/model/up0501.vtk");
    srepReader->Update();
    vtkSmartPointer<vtkPolyData> upSpokes = srepReader->GetOutput();
    vtkSmartPointer<vtkPolyData> newUpSpokes_poly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> newUpSpokes = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> spokeLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> surfacePoints = upSpokes->GetPoints();
    for(int i = 0; i < upSpokes->GetNumberOfPoints(); ++i)
    {
        double p[3];
        upSpokes->GetPoint(i, p);

        PointType hub;
        hub[0] = p[0];
        hub[1] = p[1];
        hub[2] = p[2];
        PointType transHub = tps->TransformPoint(hub);
        double newP[3];
        newP[0] = transHub[0];
        newP[1] = transHub[1];
        newP[2] = transHub[2];
        surfacePoints->SetPoint(i, newP[0], newP[1], newP[2]);
        surfacePoints->Modified();
        // int id0 = newUpSpokes->InsertNextPoint(newP);

        // double p_bdry[3];
        // upSpokes->GetPoint(i+1, p_bdry);
        // PointType bdry;
        // bdry[0] = p_bdry[0];
        // bdry[1] = p_bdry[1];
        // bdry[2] = p_bdry[2];
        // PointType transB = tps->TransformPoint(bdry);
        // double newB[3];
        // newB[0] = transB[0];
        // newB[1] = transB[1];
        // newB[2] = transB[2];
        // int id1 = newUpSpokes->InsertNextPoint(newB);

        // vtkSmartPointer<vtkLine> arrow = vtkSmartPointer<vtkLine>::New();
        // arrow->GetPointIds()->SetId(0, id0);
        // arrow->GetPointIds()->SetId(1, id1);
        // spokeLines->InsertNextCell(arrow);
    }

    newUpSpokes_poly->SetPoints(newUpSpokes);
    newUpSpokes_poly->SetLines(spokeLines);

    vtkSmartPointer<vtkPolyDataWriter> srepWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    srepWriter->SetFileName("/tmp/Slicer-zhiy/model/test.vtk");
    srepWriter->SetInput(upSpokes);
    srepWriter->Update();

    return 1;
}