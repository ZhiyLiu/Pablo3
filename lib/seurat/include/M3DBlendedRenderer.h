#ifndef M3D_BLENDED_RENDERER_H
#define M3D_BLENDED_RENDERER_H
#include <iostream>
/*
enum M3DSurfaceStyle
{
    M3D_NONE,
    M3D_WIREFRAME,
    M3D_POINTCLOUD,
    M3D_SOLID
};
//*/
#include "M3DObjectSurfaceRenderer.h"
#include <vector>
#include "Pointlist_serverB.h"
class M3DQuadFigure;
class M3DObject;
//class ThallCode::Pointlist_server2;

class Pointlist_server2;
class Diatomgrid;
class SelectedPartialFigures;

class M3DBlendedRenderer
{
public:
	M3DBlendedRenderer(M3DObject * obj = NULL, int subLevel = 1);
	virtual ~M3DBlendedRenderer();

    Diatomgrid * convertFigure(M3DQuadFigure * fig);

    void setSurfaceStyle(M3DSurfaceStyle style) { surfaceStyle = style; }
    void setPartialSurfaceStyle(M3DSurfaceStyle style) { partialSurfaceStyle = style; }
    void setLineWidth(int width) { lineWidth = width; }

	// Functions called from render(); normally not used directly
	void renderFigure(int index);
	Pointlist_server2 * computeFigure(int index, bool ignoreVisibility,
		bool ignoreBoundary);

	// Render all surfaces.  If feedback is true, nothing is drawn and the getPList()
	// function can be used to access the generated coordinate lists.  The last
	// two arguments are only used when feedback is true.  If ignoreVisibility is false,
	// only visibile figures will be used.  If ignoreBoundary is false, any existing
	// boundary displacements will be added to the resulting coordinates.
	bool render();

	void setFastInterp(int val)
	{
		std::cout << val << std::endl;
		fastInterp = val;
		std::cout << fastInterp << std::endl;
	}

	void partial_render(SelectedPartialFigures * figure_list, double cull_distance,
		bool complete_figures, int lvl);

	// Functions to access the point lists generated by render(true)
	//ThallCode::Pointlist_server2 * getPList(int figureId);

	// Function to place dots on the rendered surface.  The list is taken from
	// the figures to be rendered.  If size is greater than zero, a 6-pointed
	// jack is drawn with diameter 2*size, instead of coloring a single pixel.
	// if distinguishMarked is true, any marked landmarks will be drawn as
	// octahedrons to distinguish them, with line widths of markedWidth.
	// If the OpenGL color is not given, it will be set to contrast with the
	// figural color.

//	void drawDots(double size = 0.0, float defaultWidth = 1.0f,
//		bool distinguishMarked = false, float markedWidth = 3.0f,
//		const float * color = NULL);

	void setInterruptCallback(void (* cbFunction)(void)) { check = cbFunction; }
    //bool render();
    void interrupt() { halt = true; }

#ifndef PRODUCTION_VERSION
	// Functions for display of the blending region
    void setCurrentScale(double scale) {
		currentScale = scale;
	}
	void showBlendingRegion(bool onOff) { showRegion = onOff; }
#endif

    void init(M3DObject * obj, int subLevel = 1);
	//virtual ~M3DBlendedRenderer();

	void drawBlendedObject();

private:

	// Function to place dots on the rendered surface.  The list must be an array
	// of x,y,z coordinates of length 3*numSpots.  If size is greater than zero,
	// a 6-pointed jack is drawn with diameter 2*size, instead of coloring a single
	// pixel.  If the OpenGL color is not given, it must be set externally before
	// calling this function.  NumBigSpot may be used to specify a dot to be
	// rendered differently from the others.  It's value must refer to the
	// X-coordinate of the dot.
//	void drawDots(double * spot_list, int numSpots, double size = 0.0,
//		float defaultWidth = 1.0f, float markedWidth = 3.0f,
//		const float * color = NULL, int numBigSpot = -1);
	int lineWidth;
	int fastInterp;

	M3DObject * object;
	M3DObject * subObject;
	int level;
	int numTrees;

    M3DSurfaceStyle surfaceStyle;
    M3DSurfaceStyle partialSurfaceStyle;

//	ThallCode::Pointlist_server2 ** list;
//	int list_length;

//    Matrix4D transformMatrix;

	void (* check)(void);
//    void setSurfaceStyle(M3DSurfaceStyle style) { surfaceStyle = style; }

#ifndef PRODUCTION_VERSION
	double currentScale;	// Needed to draw blending region
	bool showRegion;		// Control for display of blending region
#endif

	static bool halt;
    //static std::vector<vertexNormal> * vertexNormalLists;
	//static std::vector<Intersection *> blendRegions;	// In breadth-first order
	static std::vector<Pointlist_server2 *> tileLists; // Indexed by figure id
//    static float * blendExtents;	// Indexed by figure id
//    static float * blendAmounts;	// Indexed by figure id

	Pointlist_serverB * blendedPList;
};

#endif

