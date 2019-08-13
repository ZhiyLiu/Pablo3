#include "M3DSpoke.h"
#include "M3DFigure.h"
#include <vnl/vnl_double_3x1.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vtkSmartPointer.h>
#include <vtksrepinterpolatecrestspokesquartic.h>
#include <vtksrepinterpolatemedialspokeshermite.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>

class M3DObject;
class M3DInterpolater
{
    public:
	
	M3DInterpolater();
	M3DInterpolater(M3DFigure *figure);
	// Method that returns a single spoke at (u,v)
	//static M3DSpoke* interpolateSpoke(M3DFigure *figure, double u, double v, int side);
	
	// Methods that return spokes given a context - usually, you'll want to use these
	// With figureId, returns all interpolated spokes for given figureId
	// Adding atomId returns all spokes around specified atom
	// Adding spokeId returns all spokes around specified spoke
	
	//static vector <M3DSpoke*> getRelevantSpokes(M3DObject* targetObject, int level, int figureId);
	//static vector <M3DSpoke*> getRelevantSpokes(M3DObject* targetObject, int level, int figureId, int atomId);
	//static vector <M3DSpoke*> getRelevantSpokes(M3DObject* targetObject, int level, int figureId, int atomId, int spokeId);
	
	// Function to draw interpolated boundary surface
	void displayInterpolatedBoundary(M3DFigure *figure, int level);
	M3DSpoke interpolateSpoke(M3DFigure* figure, double u, double v, int side);
	M3DSpoke interpolateSpoke_new( M3DFigure *figure, double u, double v, int side );
	M3DSpoke interpolateSpoke_mean(M3DFigure* figure, double u, double v, int side);
	M3DSpoke interpolateSpoke_hermite(M3DFigure* figure, double u, double v, int side);

    // interpolate 3D s-rep by new method developed by Steve and Zhiyuan in 2018
	M3DSpoke interpolate(M3DFigure *figure, double u, double v, int side);

    // interpolate 2D s-rep main entry, level is interpolate level(the power), if level == 3, then step size is (1/2)^3
    void interpolate2D(M3DFigure *figure, int level,int side, std::vector<M3DSpoke>& output);
    // interpolate recursively
    void interpolate2D(double r1, double r2, vnl_double_3x1& U1, vnl_double_3x1& U2, vnl_double_3x1& x1, vnl_double_3x1& x2, vnl_double_3x1& dx1, vnl_double_3x1& dx2, int level, double lamda, std::vector<M3DSpoke>& output);
    // interpolate skeleton
    void interpolate2D(vnl_double_3x1& x1, vnl_double_3x1& x2, double d, vnl_double_3x1& dx1, vnl_double_3x1& dx2, vnl_double_3x1& x_interp);
    // get one spoke from interpolation for 2d srep
    // input: a position want to interp
    // input: side
    // input: figure
    // output: the spoke
    void getSpokeFromInterpolation(M3DFigure* figure, int side, double x, M3DSpoke& output);
    
    // return spoke length at interpolation position
    // Meanwhile, output interpolatedU as unit dir of the spoke
    double interpolateS(M3DFigure *figure, M3DPrimitive* Sp11, M3DPrimitive* Sp12, M3DPrimitive* Sp21, M3DPrimitive* Sp22, double u, double v, int side, double lamda, Vector3D& interpolatedU);
    
    // output position of interpolated spoke
    void interpolateX(M3DFigure *figure,  int uBase, int vBase, double u, double v, int side, vnl_double_3x1& output);
    private:
	
	M3DFigure *_figure;

	// Helper functions for Hermite interpolation
    
	double h1(double s) { return 2*(s * s * s) - 3*(s * s) + 1; }
	double h2(double s) { return -2*(s * s * s) + 3*(s * s); }
	double h3(double s) { return (s * s * s) - 2*(s * s) + s; }
	double h4(double s) { return (s * s * s) - (s * s); }
	
	double h1p(double s) { return 6*(s * s) - 6*(s); }
	double h2p(double s) { return -6*(s * s) + 6*(s); }
	double h3p(double s) { return 3*(s * s) - 4*(s) + 1; }
	double h4p(double s) { return 3*(s * s) - 2*(s); }

	// Functions involved in the quaternion interpolation

	void qlog(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> &result);
	void qpow(vnl_matrix_fixed<double,4,1> q, double t, vnl_matrix_fixed<double,4,1> &result);
	void qexp(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> &result);
	void qmul(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> p, vnl_matrix_fixed<double,4,1> &result);
	void qcon(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> &result);

	vnl_matrix_fixed<double,4,1> qlog(vnl_matrix_fixed<double,4,1> q);
	vnl_matrix_fixed<double,4,1> qpow(vnl_matrix_fixed<double,4,1> q, double p);
	vnl_matrix_fixed<double,4,1> qexp(vnl_matrix_fixed<double,4,1> q);
	vnl_matrix_fixed<double,4,1> qmul(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> p);
	vnl_matrix_fixed<double,4,1> qcon(vnl_matrix_fixed<double,4,1> q);


	void slerp(vnl_double_3x1 U1, vnl_double_3x1 U2, double u, vnl_double_3x1 &result);
	vnl_matrix_fixed<double,4,1> qslerp(vnl_matrix_fixed<double,4,1> q1, vnl_matrix_fixed<double,4,1> q2, double h);
	void bislerp(vnl_double_3x1 U11, vnl_double_3x1 U12, vnl_double_3x1 U21, vnl_double_3x1 U22, double u, double v, vnl_double_3x1 &result);
	void squad(vnl_matrix_fixed<double,4,1> q0, vnl_matrix_fixed<double,4,1> q1, vnl_matrix_fixed<double,4,1> q2, vnl_matrix_fixed<double,4,1> q3, vnl_matrix_fixed<double,4,1> &result);
	vnl_double_3x1 squad(vnl_double_3x1 u0, vnl_double_3x1 u1, vnl_double_3x1 u2, vnl_double_3x1 u3, double h);
	vnl_double_3x1 dsquad(vnl_double_3x1 u0, vnl_double_3x1 u1, vnl_double_3x1 u2, vnl_double_3x1 u3, double h);
	vnl_matrix_fixed<double,4,1> g(vnl_matrix_fixed<double,4,1> qi, vnl_matrix_fixed<double,4,1> qip, vnl_matrix_fixed<double,4,1> si, vnl_matrix_fixed<double,4,1> sip, double h);
	vnl_matrix_fixed<double,4,1> dg(vnl_matrix_fixed<double,4,1> qi, vnl_matrix_fixed<double,4,1> qip, vnl_matrix_fixed<double,4,1> si, vnl_matrix_fixed<double,4,1> sip, double h);


	// Helper functions for converting Pablo types to VNL and VNL to VNL
	void convertVector3DToVNL ( Vector3D vec, vnl_double_3x1 &out );
    void convertVNLToVector3D(vnl_double_3x1 input, Vector3D& out);
	vnl_double_3x1 convertVNL4to3(vnl_matrix_fixed<double,4,1> q);
	vnl_matrix_fixed<double,4,1> convertVNL3to4(vnl_double_3x1 q);

    double getR(M3DPrimitive* p, int side);
    void getU(M3DPrimitive* p, int side, vnl_double_3x1 &u);
    
    void setR(M3DPrimitive* p, int side, double r);
    void setU(M3DPrimitive* p, int side, vnl_double_3x1 u);
    double computeMiddleS(vnl_double_3x1 startU, vnl_double_3x1 endU, vnl_double_3x1 startS,vnl_double_3x1 endS , double d, vnl_double_3x1& uMiddle);
    void compute2ndDerivative(vnl_double_3x1 startU, vnl_double_3x1 endU, vnl_double_3x1 targetU, double d, vnl_double_3x1& out);
    void printInorder(int level, std::vector<M3DSpoke>& input, std::vector<M3DSpoke>& output);
    // Description: interpolate crest spokes
    // Input: interpolation level
    // Output: crestSpokes
//    void GetCrestSpokes(vtkSmartPointer<vtkSRep> srepfig, int level, vector<M3DSpoke*>& spokes, bool istube, vtkIdType atomId, int spokeId = -1);

    // Description: generate vtk srep data structure
    // Input: quad figure
    // Output: vtk srep
    // Return: vtk id
//    vtkIdType GetVtkSrepFig(M3DFigure* figure, vtkSmartPointer<vtkSRep>& srepfig, int atomId=-1);

    };
