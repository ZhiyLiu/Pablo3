/* The purpose of this class is to apply tps on new srep representation
 *
 * Zhiyuan Liu
 * 2018.05
 * Input: folder path that host new srep files (temp.xml, up.vtp, down.vtp, crest.vtp)
 * Input: tps file path that from computePairWiseTPS
 * Input: surface mesh of ellipsoid (final shape of forward MCF)
 * Input: output prefix
 * Output: new srep files  (new_header.xml, new_up.vtp, new_down.vtp, new_crest.vtp)
 */

// Current version of new srep header format like this:
// <?xml version="1.0" ?>
// <s-rep>
// 	<nRows>5</nRows>
// 	<nCols>9</nCols>
// 	<meshType>Quad</meshType>
// 	<color>
// 		<red>0</red>
// 		<green>0.5</green>
// 		<blue>0</blue>
// 	</color>
// 	<isMean>False</isMean>
// 	<meanStatPath />
// 	<upSpoke>/playpen/workspace/demo_multiobject/802785/vtkMCF_thalamus/testNewsrep/new_down.vtp</upSpoke>
// 	<downSpoke>/playpen/workspace/demo_multiobject/802785/vtkMCF_thalamus/testNewsrep/new_down.vtp</downSpoke>
// 	<crestSpoke>/playpen/workspace/demo_multiobject/802785/vtkMCF_thalamus/testNewsrep/new_down.vtp</crestSpoke>
// </s-rep>

#include <iostream>
#include <fstream>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkSphereSource.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLDataParser.h>
#include "thinplatesplinepdmtosrep.h"
//#include "itkTransformFactoryBase.h"
//#include "itkTransformFactory.h"
//#include "itkTransformFileReader.h"
#include "itkThinPlateSplineExtended.h"


int transformNOutput(itkThinPlateSplineExtended::Pointer tps, vtkPolyData* spokes, const std::string& outputFileName)
{

    typedef double CoordinateRepType;
    typedef itkThinPlateSplineExtended TransformType;
	typedef itk::Point< CoordinateRepType, 3 > PointType;
	typedef std::vector< PointType > PointArrayType;
	typedef TransformType::PointSetType PointSetType;
	typedef PointSetType::Pointer PointSetPointer;
	typedef PointSetType::PointIdentifier PointIdType;

    vtkPointData* ptData = spokes->GetPointData();

    if(!ptData)
    {
        return -1;
    }

    vtkSmartPointer<vtkDoubleArray> arr_length = vtkDoubleArray::SafeDownCast(ptData->GetArray("spokeLength"));
    vtkSmartPointer<vtkDoubleArray> arr_dirs = vtkDoubleArray::SafeDownCast(ptData->GetArray("spokeDirection"));

	thinplatesplinepdmtosrep obj;
    vtkPoints* newPoints = vtkPoints::New();
    newPoints->SetDataTypeToDouble();
    for(int i = 0; i < spokes->GetNumberOfPoints(); ++i)
    {
        double p[3];
        spokes->GetPoint(i,p);

        // transform medial point by tps
        PointType hub;
        hub[0] = p[0];
        hub[1] = p[1];
        hub[2] = p[2];
        PointType transHub = tps->TransformPoint(hub);
        double newP[3];
        newP[0] = transHub[0];
        newP[1] = transHub[1];
        newP[2] = transHub[2];
        newPoints->InsertNextPoint(newP);
        
        // transform implied boundary point by tps
        PointType bdryPt;
        int baseIdx = i * 3;
        double length = arr_length->GetValue(i);
        double dirX = arr_dirs->GetValue(baseIdx);
        double dirY = arr_dirs->GetValue(baseIdx+1);
        double dirZ = arr_dirs->GetValue(baseIdx+2);
        bdryPt[0] = p[0] + dirX * length;
        bdryPt[1] = p[1] + dirY * length;
        bdryPt[2] = p[2] + dirZ * length;
        PointType transBdryPt = tps->TransformPoint(bdryPt);

        double newLength = obj.calculateSpokeLength(transHub, transBdryPt);
        arr_length->SetValue(i, newLength);

        Vector3D newDir = obj.calculateSpokeDirection(transHub, transBdryPt);
        arr_dirs->SetValue(baseIdx, newDir.getX());
        arr_dirs->SetValue(baseIdx+1, newDir.getY());
        arr_dirs->SetValue(baseIdx+2, newDir.getZ());

    }
    spokes->SetPoints(newPoints);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetDataModeToAscii();
    writer->SetFileName(outputFileName.c_str());
    writer->SetInput(spokes);
    writer->Update();

    return 1;
}

// this function assume the xml file has 3 levels( root is level 1) at most
void writeHeader(vtkXMLDataElement* root, const std::string& outputPrefix)
{
    std::ofstream headerFile;
    std::string fileName = outputPrefix+"new_header.xml";
    headerFile.open(fileName.c_str());
    headerFile << "<?xml version=\"1.0\" ?>\n";
    headerFile << "<" << root->GetName() << ">\n";
    for(int i = 0; i < root->GetNumberOfNestedElements(); ++i)
    {
        vtkXMLDataElement* component = root->GetNestedElement(i);
        std::string name = component->GetName();
        std::string value = component->GetCharacterData();
        if(component->GetNumberOfNestedElements() == 0 && value.empty())
        {
            // empty node
            headerFile << "\t<" << name << " \/>\n";
        }
        else if(component->GetNumberOfNestedElements() == 0)
        {
            // attribute
            headerFile << "\t<" << name << ">" << value << "<\/" << name << ">\n";
        }
        else
        {
            // have children
            headerFile << "\t<" << name << ">\n";
            for(int j = 0; j < component->GetNumberOfNestedElements(); ++j)
            {
                vtkXMLDataElement* sub_component = component->GetNestedElement(j);
                std::string myName = sub_component->GetName();
                std::string myValue = sub_component->GetCharacterData();
                headerFile << "\t\t<" << myName << ">" << myValue << "<\/" << myName << ">\n";
            }
            headerFile << "\t<\/" << name << ">\n";
        }
    }
    headerFile << "<\/" << root->GetName() << ">\n";
    headerFile.close();
}
int main(int argc, char** argv)
{
	if( argc != 5 ) {
		std::cerr << "Usage: "<< std::endl;
		std::cerr << argv[0];
		std::cerr << " <Folder of new s-rep, tpsFileName, sourceLandMark and outputPrefix>";
		std::cerr << std::endl;
		return -1;
	}
	std::string templateModelName=argv[1];
	std::string transformFileName=argv[2];
	std::string sourceLandMarkFileName=argv[3];
	std::string outputPrefix=argv[4];

    std::string outputUpFileName = outputPrefix + "new_up.vtp";
    std::string outputDownFileName = outputPrefix + "new_down.vtp";
    std::string outputCrestFileName = outputPrefix + "new_crest.vtp";

    // 1. read in new srep format
    // Get all data from the file up spokes
    vtkSmartPointer<vtkXMLDataParser> headerParser = vtkSmartPointer<vtkXMLDataParser>::New();
    std::string headerFileName = templateModelName + "temp.xml";
    headerParser->SetFileName(headerFileName.c_str());

    if(headerParser->Parse() != 1)
    {
        std::cout << "Failed to parse header file." << std::endl;
        return -1;
    }
    vtkSmartPointer<vtkXMLDataElement> root = headerParser->GetRootElement();
    int numOfElements = root->GetNumberOfNestedElements();
    std::vector<std::string> vtpFiles;
    std::vector<std::string> outputFiles;
    outputFiles.push_back(outputUpFileName);
    outputFiles.push_back(outputDownFileName);
    outputFiles.push_back(outputCrestFileName);
    for(int i = 0; i < numOfElements; ++i)
    {
        vtkSmartPointer<vtkXMLDataElement> component = root->GetNestedElement(i);
        std::string name = component->GetName();
        std::size_t found = name.find("Spoke");
        int fileID = 0;
        if(found != std::string::npos)
        {
            fileID++;
            vtpFiles.push_back(component->GetCharacterData());
            component->SetCharacterData(outputFiles[fileID].c_str(), outputFiles[fileID].length());
        }
    }
    
    // assume the 1st is up, 2nd is down, 3rd file is crestSpokes
    if(vtpFiles.size() < 3)
    {
        std::cout << "vtp files should be up.vtp, down.vtp, crest.vtp" << std::endl;
        return -1;
    }
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(vtpFiles[0].c_str());
    reader->Update();
    vtkPolyData* upSpokes = reader->GetOutput();

    // Get down spokes
    vtkSmartPointer<vtkXMLPolyDataReader> reader2 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader2->SetFileName(vtpFiles[1].c_str());
    reader2->Update();
    vtkPolyData* downSpokes = reader2->GetOutput();

    // Get crest
    vtkSmartPointer<vtkXMLPolyDataReader> reader3 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
//    std::string crestFileName = templateModelName + "crest.vtp";
    reader3->SetFileName(vtpFiles[2].c_str());
    reader3->Update();
    vtkPolyData* crestSpokes = reader3->GetOutput();

    typedef double CoordinateRepType;
    typedef itkThinPlateSplineExtended TransformType;
	typedef itk::Point< CoordinateRepType, 3 > PointType;
	typedef std::vector< PointType > PointArrayType;
	typedef TransformType::PointSetType PointSetType;
	typedef PointSetType::Pointer PointSetPointer;
	typedef PointSetType::PointIdentifier PointIdType;

    // 2. prepare for transformation matrix
	ifstream inFile;
	inFile.open(transformFileName.c_str());
	if(!inFile) {
		std::cerr << "Unable to open the file: " << transformFileName << std::endl;
		return EXIT_FAILURE;
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
	reader_source->SetFileName(sourceLandMarkFileName.c_str());
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

    // 3. transform and output
    // iterate up spokes
    transformNOutput(tps, upSpokes, outputUpFileName);

    // process down spokes
    transformNOutput(tps, downSpokes, outputDownFileName);

    // process crest spokes
    transformNOutput(tps, crestSpokes, outputCrestFileName);

    // output header file
    writeHeader(root, outputPrefix);
    return 1;
}