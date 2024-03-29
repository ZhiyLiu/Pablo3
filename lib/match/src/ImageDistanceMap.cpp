//#include <sys/timeb.h>
//#include <unistd.h>

#include <iostream>
#include <math.h>
#include "Danielsson.h"
#include "DistanceMap3D.h"
#include "ImageDistanceMap.h"
#include "ImageResample3D.h"
#include "snapshot.h"
#include "Tuning.h"



//#define DEBUG


extern int globalVerbosity;

using namespace std;


ImageDistanceMap::ImageDistanceMap()
{
	distIm = NULL;
	//xiaojie
	gradDistIm = NULL;
}


ImageDistanceMap::ImageDistanceMap(Image3D * im)
{
	image = im;

	distIm = NULL;
	//xiaojie
	gradDistIm = NULL;

	initialize();
}


ImageDistanceMap::~ImageDistanceMap()
{
	if (distIm != NULL)
		delete [] distIm;
	if (image != NULL)
		delete image;
}


const char * msg_1 = "Error: ImageDistanceMap object not initialized";


void ImageDistanceMap::initialize(Image3D * im)
{
	if (im == NULL) {
		cout << msg_1 << endl;
		return;
	}
	if (distIm != NULL) {
		delete [] distIm;
		distIm = NULL;
	}
	image = im;
	initialize();
}


void ImageDistanceMap::initialize()
{
#ifdef BINARY
	if (image->getAbsXSpacing() != image->getAbsYSpacing()
		|| image->getAbsXSpacing() != image->getAbsZSpacing())
	{
		// Deal with anisotropy by just resampling image to isotropic.
		// Easy, but expensive in time.
		if (globalVerbosity > -1) {
			cout << "Resampling image to cubic voxels" << endl;
		}
		ImageResample3D resampler;
		resampler.isotropicSample(*image, tuningWt(BpVoxelSpacing), 
			tuningWt(AttractionMask), tuningWt(RepulsionMask));
	}
#endif

	mapSize[0] = image->getXDim();
	mapSize[1] = image->getYDim();
	mapSize[2] = image->getZDim();
	size = mapSize[0] * mapSize[1] * mapSize[2];

	// Or deal with anisotropy by just scaling everything properly
	double shortest_spacing = min(image->getAbsXSpacing(), min(image->getAbsYSpacing(),
		image->getAbsZSpacing()));
	scale = 1.0/shortest_spacing;

	// This puts all distances into units of 'shortest-side-voxels'
	mapSpacing[0] = image->getAbsXSpacing() * scale;
	mapSpacing[1] = image->getAbsYSpacing() * scale;
	mapSpacing[2] = image->getAbsZSpacing() * scale;

	const double extent[] = {
		mapSize[0]*mapSpacing[0]/scale,
		mapSize[1]*mapSpacing[1]/scale,
		mapSize[2]*mapSpacing[2]/scale
	};

	// The scaling from model to world coordinates is the maximum extent
	double modelToWorldScale = extent[0];
	if (modelToWorldScale < extent[1])
		modelToWorldScale = extent[1];
	if (modelToWorldScale < extent[2])
		modelToWorldScale = extent[2];

	// Scaling from model to image coordinates
	modelToImageScale[0] = modelToWorldScale / mapSpacing[0] * scale;
	modelToImageScale[1] = modelToWorldScale / mapSpacing[1] * scale;
	modelToImageScale[2] = modelToWorldScale / mapSpacing[2] * scale;

	if (globalVerbosity > 0) {
		cout << "Mapsize: " << mapSize[0] << ", " << mapSize[1] << ", " 
			<< mapSize[2] << endl;
		cout << "Mapspacing: " << mapSpacing[0] << ", " << mapSpacing[1] << ", " 
			<< mapSpacing[2] << endl;
	}
}

// return count of 6-connected voxel neighbors that are off (non-positive).
// The (x,y,z) coordinate is in a *raster* order that ignores image coordinate
// flipping and is not an array index only to avoud wrapping around the image edge.
inline int ImageDistanceMap::zeroNeighbors(int x, int y, int z)
{
	int count = 0;

	int xdim = mapSize[0];
	int ydim = mapSize[1];
	int zdim = mapSize[2];
	int xydim = xdim*ydim;
	int pixIndex = x + xdim*y + xydim*z;
	GreyValue *vox = image->getVoxels();

	// ignore neighbors that wrap around the image's edge
	if (x           && vox[pixIndex-    1] <= 0) count++;
	if (y           && vox[pixIndex- xdim] <= 0) count++;
	if (z           && vox[pixIndex-xydim] <= 0) count++;
	if (x < xdim-1  && vox[pixIndex+    1] <= 0) count++;
	if (y < ydim-1  && vox[pixIndex+ xdim] <= 0) count++;
	if (z < zdim-1  && vox[pixIndex+xydim] <= 0) count++;

	return count;
}


// Compute the distance map for the binary image.

bool ImageDistanceMap::createMap()
{
	int x, y, z;
	int i;

	int boundary_count = 0;

	if (image == NULL || size == 0) {
		cout << msg_1 << endl;
		return false;
	}

	unsigned char * bdryIm = new unsigned char[size];
	// use bitset instead?

	distIm = new short[size];

	// should use memset() here
	for (i = 0; i < size; i++) {
		bdryIm[i] = 0;
		distIm[i] = 0;
	}

	// Convert the binary/blurred solid object image into a
	// boundary (shell) object image. The image is blurred if
	// resampling was done by initialize(), else it's derived from
	// a plane of a bitStack image or just a binary image (which has
	// only one plane).

	// Image flip is ignored throughout distance map calc's:
	// avoid Image3D voxel accessors by accessing the voxels directly.
	// This is ~6x faster and bdryIm is 2x smaller.
	// --GST 20080111
#ifdef UNFLIPPED
	cout
		<< "Warning: distance map code is untested if UNFLIPPED is defined at compile time."
		<< endl;
#endif

	GreyValue *binImVox = image->getVoxels();
	int pixIndex = 0;
	for (z = 0; z < mapSize[2]; z++) {
		for (y = 0; y < mapSize[1]; y++) {
			for (x = 0; x < mapSize[0]; x++, pixIndex++) {

				// If this voxel is on and has no neighbors, it's a boundary.
				// Threshold at 0.5: works for blurred OR binary image.
				// NOTE: zeroNeighbors assumes array (x,y,z) -- no flipping
				if (binImVox[pixIndex] > 0.5 && zeroNeighbors(x, y, z)) {
					bdryIm[pixIndex] = (unsigned char) 255;
					boundary_count++;
				}
			}
		}
	}

#ifdef BINARY
	if ((int) tuningWt(BpDebug1) == 1) {
		cout << "Snapshot bdryIm.raw3" << endl;
		snapshot("bdryIm.raw3", bdryIm, mapSize[0], mapSize[1], mapSize[2],
			0, 255, mapSpacing[0], mapSpacing[1], mapSpacing[2]);
	}
#endif

	if (globalVerbosity > 0)
		cout << "Boundary points found: " << boundary_count << endl;

	bool ret = createRawSignedImageDistanceMaps((char *) bdryIm);
	delete [] bdryIm;
	return ret;
}

//#define ROHIT_DEBUG
#ifdef ROHIT_DEBUG
#include <fstream>
#endif


bool ImageDistanceMap::createRawSignedImageDistanceMaps(char * bdryIm)
{
	short * odx;
	short * ody;
	short * odz;
	int i;
	float dist;

#ifdef DEBUG
	float dmin = 1e5;
	float dmax = -1e5;
#endif

	if (globalVerbosity > -1) {
		cout << "Calculating Danielsson Distance Map (DDM)" << endl;
	}
#ifdef ROHIT_DEBUG
	ofstream fbdry("bdry.dat", ios::trunc|ios::out);
	for (i = 0; i < size; i++)
	{
		fbdry << (unsigned short)bdryIm[i] << endl;
	}
#endif

	int retCode;	// failure code for edt3ddan; 0=OK, else=outOfMemory
	retCode = edt3ddan(bdryIm, mapSize[0], mapSize[1], mapSize[2],
#ifdef BINARY
		mapSpacing[0], mapSpacing[1], mapSpacing[2],
#endif
		0, &odx, &ody, &odz);
	if (retCode)
	{
		cout << "Error: Insufficient memory for distance map"
			<< ", err code=" << retCode
			<< ", mapsize=" << mapSize[0]*mapSize[1]*mapSize[2]*2 << " *3*shorts"
			<< endl;
		return false;
	}

#ifdef DEBUG
	cout << "createRawSignedDistanceBpointMaps(): size = " << size << endl;
#endif


	// Deal with anisotropy in units of 'shortest-side-voxels'
	double x_spacing_sqare = mapSpacing[0]*mapSpacing[0];
	double y_spacing_sqare = mapSpacing[1]*mapSpacing[1];
	double z_spacing_sqare = mapSpacing[2]*mapSpacing[2];

#ifdef ROHIT_DEBUG
	cout << mapSize[0] << ' ' << mapSize[1] << ' ' << mapSize[2] << endl;

	ofstream fmap("dmap.dat", ios::trunc|ios::out);
	ofstream fimg("img.dat", ios::trunc|ios::out);
#endif

	GreyValue *binaryImage = image->getVoxels();

	//#define DEBUG_DDM
#ifdef DEBUG_DDM
	short testx[] = {0, 1, 1, 1, 2, 0, 1, 1, 1, 2};
	short testy[] = {0, 0, 1, 1, 0, 0, 0, 1, 1, 0};
	short testz[] = {0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
	short binIm[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
#endif

	GreyValue stackMask = (image->getIsImageStacked() ? image->getStackedMask() : 0xFFFF);
	for (i = 0; i < size; i++)
	{
#ifdef DEBUG
		if( i % 10000 == 0 ) {
			printf("\b\b\b\b\b\b %3d%% ", int(double(i)*100./double(size)) );
			fflush(stdout);
		}
#endif
		// Okay for isotropic, ie, if we resample
		//dist = (float) (sqrt(((double) odx[i]*odx[i] + ody[i]*ody[i] + odz[i]*odz[i])));
		// Req'd for anisotropic spacing, can't hurt - with isotropic voxels, the scale is just 1:1

		//dist = (float) (sqrt(((double) odx[i]*odx[i]*x_spacing_sqare +
		//							   ody[i]*ody[i]*y_spacing_sqare +
		//							   odz[i]*odz[i]*z_spacing_sqare  )));

		// GST: test moving the boundary from the voxel centers to the voxel edges:
		//  expand each internal distance by a half-voxel and contract the external distances.
		//  - round distance to nearest fixed-point integer, then add
		//      .5 voxel if inside the image 
		//     -.5 voxel if outside.
		//  - the goal is a +1/2 voxel distance for the voxel bordering the "voxel
		//    crack" and -1/2 voxel distance just inside the crack, with a zero halfway
		//    between the centers (at the crack).
		//  - routines that use the distance map may need to take fabs(mapValue) now.
		//  - we really want to store dist^2, but there are not enough bits in a short int.
		//    800 is the max dist we expect, so 800*800 = 640000 -- bigger than 16 bits
		//  - distances inside the object are negative; outide are positive.
#ifdef DEBUG_DDM
		if (i >= sizeof(testx)/sizeof(*testx)) exit(0);
		float dist2 =	// distance squared
			testx[i]*testx[i] +
			testy[i]*testy[i] +
			testz[i]*testz[i];
#else
		float dist2 =	// distance squared
			odx[i]*odx[i]*x_spacing_sqare +
			ody[i]*ody[i]*y_spacing_sqare +
			odz[i]*odz[i]*z_spacing_sqare;
#endif
		dist = sqrt(dist2);

#ifdef DEBUG
		if (dist > dmax)
			dmax = dist;
		if (dist < dmin)
			dmin = dist;
#endif

#ifdef BINARY
		// Store distance (voxels)
		if ((int) tuningWt(BpDDMCenter) == 1)

		{			
			// Round to nearest integer voxel distance.
			// Side effect: forces most corners of 4-connected neighbors and 8-connected
			//   neighbors to a distance of 1, eg, sqrt(1+1+0) -> 1.4 -> 1
			distIm[i] = (short) (dist + 0.5);
		}
		else {
			// We believe that this has the effect of a zero crossing halfway between the
			// inside and outside voxels.
			//distIm[i] = (short) (8.0*dist + 0.5) + (binaryImage[i] ? 4 : -4);
			//distIm[i] = (short) (8.0*dist + 0.5);	// round to nearest signed 12.3 integer voxel distance
#ifdef DEBUG_DDM
			distIm[i] = floor(8.0*dist + 0.5);	// round to nearest signed 12.3 integer voxel distance
			if (binIm[i]) distIm[i] = -distIm[i];
			distIm[i] = distIm[i] - 4.0;
			cout << "DEBUG_DDM: (x,y,z,binIm)=" << testx[i] << " " << testy[i] << " " << testz[i]
			<< " " << binIm[i] << ") " << "distIm=" << distIm[i] << endl;
#else
			distIm[i] = (short) floor(8.0*dist + 0.5);	// round to nearest signed 12.3 integer voxel distance
			if (binaryImage[i] & stackMask)  {
				distIm[i] = -distIm[i];
			}
			distIm[i] = distIm[i] - 4;
#endif
		}
#endif	/* BINARY */
#ifdef ROHIT_DEBUG
		fmap << (short) distIm[i] << endl;
		fimg << (short) binaryImage[i] << endl;
#endif
	}
	cout << '\n';

	free(odx);
	free(ody);
	free(odz);

#ifdef BINARY
	if ((int) tuningWt(BpDebug2) == 1)
	{	
		cout << "Snapshot ddm.raw3" << endl;

		float dmin = (float) (distIm[0]), dmax = (float)(distIm[0]);
		float * fltMap = new float[size];
		dmin = distIm[0];
		dmax = dmin;
		for (int k = 0; k < size; k++) {
			fltMap[k] = (float) (distIm[k]);
			if (fltMap[k] < dmin)
				dmin = fltMap[k];
			else if (fltMap[k] > dmax)
				dmax = fltMap[k];
		}

		snapshot("ddm.raw3", fltMap, mapSize[0], mapSize[1], mapSize[2],
			(float)dmin, (float)dmax, mapSpacing[0], mapSpacing[1], mapSpacing[2]);

		delete fltMap;

		cout << "Distance map intensity range = " << dmin << " - " << dmax << endl;
	}
#else
#ifdef DEBUG
	cout << "Distance map intensity range = " << dmin << " - " << dmax << endl;
#endif
#endif
	return true;
}


/*	Assuming that all distances outside are positive, this function
extrapolates the distance map using city block distances.  This
is inaccurate, but monotonically increasing, so it is appropriate
for optimization.
*/
float ImageDistanceMap::getDistance(int i, int j, int k)
{
	int inc	= 0;
	if (i < 0) {
		inc	-= i;
		i = 0;
	}
	else if (i >= mapSize[0]) {
		inc	+= (i - mapSize[0] + 1);
		i = mapSize[0] - 1;
	}
	if (j < 0) {
		inc	-= j;
		j = 0;
	}
	else if (j >= mapSize[1]) {
		inc	+= (j - mapSize[1] + 1);
		j = mapSize[1] - 1;
	}
	if (k < 0) {
		inc	-= k;
		k = 0;
	}
	else if (k >= mapSize[2]) {
		inc	+= (k - mapSize[2] + 1);
		k = mapSize[2] - 1;
	}

	double ret = (int) distIm[i + mapSize[0]*(j + k*mapSize[1])];
	
	// xiaojie: earlier it was ret/8, now it is ret/100
	//return ret/8 + inc;
	return ret/100 + inc;	// this division is done because distances are stored as unsigned shorts 
							// which are 100 times the actual distances
}


float ImageDistanceMap::bilerp_in_z(float x, float y, int z)
{
	int ix, iy;
	float   dx, dy;

	ix = (int) x;
	iy = (int) y;
	float a = getDistance(ix, iy, z);
	float b = getDistance((ix + 1), iy, z);
	float c = getDistance(ix, (iy + 1), z);
	float d = getDistance((ix + 1), (iy + 1), z);
	dx = x - (float) ix;
	dy = y - (float) iy;

	float i1 = (1.0 - dx)*a + dx*b;
	float i2 = (1.0 - dx)*c + dx*d;
	return (1.0 - dy)*i1 + dy*i2;
}


float ImageDistanceMap::trilerp(float x, float y, float z)
{
	int    ix, iy, iz;
	float   a, b;
	float   dz;

	ix = (int) x;
	iy = (int) y;
	iz = (int) z;

	a = bilerp_in_z(x, y, iz);
	b = bilerp_in_z(x, y, iz + 1);
	dz = z - (float) iz;

	return (1.0 - dz)*a + dz*b;
}


const char * msg_2 =
"Error: ImageDistanceMap object not initialized: returning infinite distance";


// This function doesn't test to verify that (x, y, z) is in the unit cube
float ImageDistanceMap::get_distance(double x, double y, double z)
{
	if (distIm == NULL) {
		if (globalVerbosity > 1)
			cout << msg_2 << endl;
		return INFINITE_DISTANCE;
	}
	Vector3D v(x, y, z);
	image->modelToImageCoordinates(v);
	// Points outside the unit cube will return monotonically
	// increasing values, by use of the city-block distance
	return trilerp(v.getX(), v.getY(), v.getZ());
}
//xiaojie
//doesn't use this function in current program-11/22/2011
Vector3D ImageDistanceMap::getGradDistance(int i, int j, int k)
{



	//TOCHECK
	//return gradDistIm[i + mapSize[0]*(j + k*mapSize[1])];

	return( Vector3D() ) ;

}


Vector3D ImageDistanceMap::bilerp_in_z_grad(float x, float y, int z)
{
	int ix, iy;
	float   dx, dy;

	ix = (int) x;
	iy = (int) y;
	Vector3D a = getGradDistance(ix, iy, z);
	Vector3D b = getGradDistance((ix + 1), iy, z);
	Vector3D c = getGradDistance(ix, (iy + 1), z);
	Vector3D d = getGradDistance((ix + 1), (iy + 1), z);
	dx = x - (float) ix;
	dy = y - (float) iy;

	Vector3D i1 = a*(1.0 - dx) + b*dx;
	Vector3D i2 = c*(1.0 - dx) + d*dx;
	return i1*(1.0 - dy) + i2*y;
}


Vector3D ImageDistanceMap::trilerp_grad(float x, float y, float z)
{
	int    ix, iy, iz;
	Vector3D   a, b;
	float   dz;

	ix = (int) x;
	iy = (int) y;
	iz = (int) z;

	a = bilerp_in_z_grad(x, y, iz);
	b = bilerp_in_z_grad(x, y, iz + 1);
	dz = z - (float) iz;

	return a*(1.0 - dz) + b*dz;
}


//const char * msg_2 =
//"Error: ImageDistanceMap object not initialized: returning infinite distance";


// This function doesn't test to verify that (x, y, z) is in the unit cube
Vector3D ImageDistanceMap::get_graddistance(double x, double y, double z)
{

	// ----------------------------------------------------------

	const Vector3D imageToModelScale(
		1.0/getModelToImageScale()[0],
		1.0/getModelToImageScale()[1],
		1.0/getModelToImageScale()[2] );

	// dp is delta voxel in unit co-ordinates.
	Vector3D dp = imageToModelScale;
	dp = dp / 2.;	// halve the step size as an initial guess
	double dpx = dp.getX();		// faster forms of dp
	double dpy = dp.getY();
	double dpz = dp.getZ();

	// The gradient is derived from the distance map, a piecewise-linear map which is tri-lerp'd by
	// getDistance() at positions slightly around (x,y,z). The second term (weighted by .0532) helps
	// decide cases in which getDistance(+dp) == getDistance(-dp).
	// FIXME: Check for bounds
	// (No need to adjust for dilation here, as we are only interested in differences)

	Vector3D grad = 0.9679 * Vector3D(
		getDistance(x+dpx,y,z) - getDistance(x-dpx,y,z),
		getDistance(x,y+dpy,z) - getDistance(x,y-dpy,z),
		getDistance(x,y,z+dpz) - getDistance(x,y,z-dpz) )
		+ 0.0532 * Vector3D(
		getDistance(x+2*dpx,y,z) - getDistance(x-2*dpx,y,z),
		getDistance(x,y+2*dpy,z) - getDistance(x,y-2*dpy,z),
		getDistance(x,y,z+2*dpz) - getDistance(x,y,z-2*dpz) );

	if(grad == Vector3D(0,0,0)){
		grad = 0.9679 * Vector3D(
			getDistance(x+5*dpx,y,z) - getDistance(x-5*dpx,y,z),
			getDistance(x,y+5*dpy,z) - getDistance(x,y-5*dpy,z),
			getDistance(x,y,z+5*dpz) - getDistance(x,y,z-5*dpz) )
			+ 0.0532 * Vector3D(
			getDistance(x+10*dpx,y,z) - getDistance(x-10*dpx,y,z),
			getDistance(x,y+10*dpy,z) - getDistance(x,y-10*dpy,z),
			getDistance(x,y,z+10*dpz) - getDistance(x,y,z-10*dpz) );
		//grad.normalize();


	}


	grad.normalize();

	return( grad ) ;

	// ----------------------------------------------------------

}

float ImageDistanceMap::getWorldDistance(double x, double y, double z)
{
	// Gets something in unit cube distances
	if (x < 0.0 || y < 0.0 || z < 0.0)
		return (float) INFINITE_DISTANCE;
	if (x >= 1.0 || y >= 1.0 || z >= 1.0)
		return (float) INFINITE_DISTANCE;
	if (distIm == NULL) {
		if (globalVerbosity > 1)
			cout << msg_2 << endl;
		return INFINITE_DISTANCE;
	}

	// Converts them to world distances
	Vector3D v(x, y, z);
	image->modelToImageCoordinates(v);

	int i = (int) v.getX();         // AGG: Why not rounding here?
	int j = (int) v.getY();
	int k = (int) v.getZ();

	return (float) image->modelToWorldDistance(getDistance(i, j, k));
}


ImageDistanceMap::operator Image3D *()
{
	// the output voxels
	GreyValue * img = new GreyValue[size];

	int y,z;
	int xdim = mapSize[0];
	int ydim = mapSize[1];
	int zdim = mapSize[2];
	int xydim = xdim * ydim;
	if (image->getXSpacing() < 0)
	{
		cout << "Error: x flip is not implemented";
	}
	else if (image->getYSpacing() < 0)
	{
		for (z=0; z < zdim; z++)
		{
			for (y=0; y < ydim; y++)
			{
				memcpy((void *) (img    + z*xydim + y*xdim),
					(void *) (distIm + z*xydim +(ydim-y-1)*xdim),
					xdim * sizeof(GreyValue));
			}
		}
	}
	else if (image->getZSpacing() < 0)
	{
		cout << "Error: z flip is not implemented";
	}
	else
	{
		// no flip - copy without changing voxel order
		memcpy((void *) img, (void *) distIm, size * sizeof(GreyValue));
	}

	Image3D * im = new Image3D;
	im->setVoxels(img, mapSize[0], mapSize[1], mapSize[2]);
	// The image will have all spacings positive
	im->setSpacingAndOrigin(mapSpacing[0]/scale, mapSpacing[1]/scale, mapSpacing[2]/scale);

	return im;
}


void ImageDistanceMap::printMinMaxDistance()
{
	float min = distIm[0];
	float max = distIm[0];
	for (int i = 0; i < size; i++) {
		if (max < distIm[i])
			max = distIm[i];
		else if (min > distIm[i])
			min = distIm[i];
	}

	cout << " distances in map: " << min << " .. " << max << endl;
}


void ImageDistanceMap::fromImage3D(Image3D * distMapImage)
{
	initialize(distMapImage);	// This deletes distIm
	distIm = new short[size];

	GreyValue * im = distMapImage->getVoxels();

	for (int i = 0; i < size; i++) {
		distIm[i] = im[i];
		//if (i%1000==0) cout << im[i] << " : " << distIm[i] << endl;
	}

	if (globalVerbosity > 0)
		printMinMaxDistance();

	return;
}
//xiaojie

void ImageDistanceMap::createGradDistMapFromGradXYZ( Image3D* gradX, Image3D* gradY, Image3D* gradZ )
{
	// delete the gradDistIm if it already exists

	if( gradDistIm != NULL ) {
		delete [] gradDistIm ;
		gradDistIm = NULL ;
	}

	gradDistIm = new Vector3D[size];

	GreyValue * imX = gradX->getVoxels();
	GreyValue * imY = gradY->getVoxels();
	GreyValue * imZ = gradZ->getVoxels();

	for (int i = 0; i < size; i++) {
		gradDistIm[i].setX(imX[i]);
		gradDistIm[i].setY(imY[i]);
		gradDistIm[i].setZ(imZ[i]);		
	}

	return;
}


void ImageDistanceMap::setgradXDistImfromImage3D(Image3D * gradXdistImage)
{
	// delete the gradDistIm pointer if it already exists


	gradDistIm = new Vector3D[size];

	GreyValue * im = gradXdistImage->getVoxels();

	for (int i = 0; i < size; i++) {
		gradDistIm[i].setX(im[i]);
		//if (i%1000==0) cout << im[i] << " : " << distIm[i] << endl;
	}

	return;
}
void ImageDistanceMap::setgradYDistImfromImage3D(Image3D * gradYdistImage)
{
	gradDistIm = new Vector3D[size];

	GreyValue * im = gradYdistImage->getVoxels();

	for (int i = 0; i < size; i++) {
		gradDistIm[i].setY(im[i]);
		//if (i%1000==0) cout << im[i] << " : " << distIm[i] << endl;
	}

	return;
}
void ImageDistanceMap::setgradZDistImfromImage3D(Image3D * gradZdistImage)
{
	gradDistIm = new Vector3D[size];

	GreyValue * im = gradZdistImage->getVoxels();

	for (int i = 0; i < size; i++) {
		gradDistIm[i].setZ(im[i]);
		//if (i%1000==0) cout << im[i] << " : " << distIm[i] << endl;
	}

	return;
}
// Xiaojie

Curvature ImageDistanceMap::getCurvature(double x, double y, double z) {

	// get the gradient of the distance

	Vector3D n = getGradDistance(x,y,z);
	//cout<<"n"<<n<<endl;

	// let(nX, nY, nZ) be the three components of the n vector

	// get an u,v basis for the tangent plane
	//TOFIX
	//n=(0,0,1)

	Vector3D u,v;

	if((fabs(n.getY())<1e-10)&&(fabs(n.getX())<1e-10)){
		u = Vector3D(1,0,0);
		v = Vector3D(0,1,0);
	}
	else{
		//xiaojie0426
		u = Vector3D(n.getY(), -n.getX(), 0);
		u.normalize();
		v = n.cross(u);
		v.normalize();
	}

	//DEBUG
	//cout<<"n"<<n<<endl;
	//cout<<"u"<<u<<endl;
	//cout<<"v"<<v<<endl;

	// For calculating the unit distance 

	const Vector3D imageToModelScale(
		1.0/getModelToImageScale()[0],
		1.0/getModelToImageScale()[1],
		1.0/getModelToImageScale()[2] );

	// dp is delta voxel in unit co-ordinates.

	Vector3D dp = imageToModelScale;
	dp = dp / 2.;	// halve the step size as an initial guess
	double dpx = dp.getX();		// faster forms of dp
	double dpy = dp.getY();
	double dpz = dp.getZ();

	// calculating the gradient of gradient of distance

	Vector3D gradGrad_X =	0.9679 * ( getGradDistance(x+dpx,y,z) - getGradDistance(x-dpx,y,z) ) 
		+	0.0532 * ( getGradDistance(x+2*dpx,y,z) - getGradDistance(x-2*dpx,y,z) )   ; 
	if(gradGrad_X==Vector3D(0,0,0)){
		gradGrad_X =	0.9679 * ( getGradDistance(x+5*dpx,y,z) - getGradDistance(x-5*dpx,y,z) ) 
			+	0.0532 * ( getGradDistance(x+10*dpx,y,z) - getGradDistance(x-10*dpx,y,z) )   ; 

	}


	Vector3D gradGrad_Y =	0.9679 * ( getGradDistance(x,y+dpy,z) - getGradDistance(x,y-dpy,z) ) 
		+	0.0532 * ( getGradDistance(x,y+2*dpy,z) - getGradDistance(x,y-2*dpy,z) )   ;
	if(gradGrad_Y==Vector3D(0,0,0)){
		gradGrad_Y =	0.9679 * ( getGradDistance(x,y+5*dpy,z) - getGradDistance(x,y-5*dpy,z) ) 
			+	0.0532 * ( getGradDistance(x,y+10*dpy,z) - getGradDistance(x,y-10*dpy,z) )   ;
	}


	Vector3D gradGrad_Z =	0.9679 * ( getGradDistance(x,y,z+dpz) - getGradDistance(x,y,z-dpz) ) 
		+	0.0532 * ( getGradDistance(x,y,z+2*dpz) - getGradDistance(x,y,z-2*dpz) )   ; 
	if(gradGrad_Z==Vector3D(0,0,0)){
		gradGrad_Z =	0.9679 * ( getGradDistance(x,y,z+5*dpz) - getGradDistance(x,y,z-5*dpz) ) 
			+	0.0532 * ( getGradDistance(x,y,z+10*dpz) - getGradDistance(x,y,z-10*dpz) )   ; 
	}


	// calculating the gradients of the gradient n = (nX, nY, nZ)

	Vector3D grad_nX = Vector3D( gradGrad_X.getX(), gradGrad_Y.getX(), gradGrad_Z.getX() ) ;
	Vector3D grad_nY = Vector3D( gradGrad_X.getY(), gradGrad_Y.getY(), gradGrad_Z.getY() ) ;
	Vector3D grad_nZ = Vector3D( gradGrad_X.getZ(), gradGrad_Y.getZ(), gradGrad_Z.getZ() ) ;

	// Curvature Tensor: the 2X2 matrix (M-II shape operator) in u,v basis

	double a11 = u * Vector3D(grad_nX*u, grad_nY*u, grad_nZ*u);
	double a12 = u * Vector3D(grad_nX*v, grad_nY*v, grad_nZ*v);
	double a22 = v * Vector3D(grad_nX*v, grad_nY*v, grad_nZ*v);

	// Construct the curvature from the tensor values	

	// calculating principal curvatures
	//TOCHECK
	//kappa1>=kapap2?

	double squaredTerm = a11*a11+4*a12*a12-2*a11*a22+a22*a22 ;

	if( squaredTerm < 0 )
		squaredTerm = 0 ;

	double kappa1 = double(0.5)*(a11+a22+sqrt(squaredTerm));
	double kappa2 = double(0.5)*(a11+a22-sqrt(squaredTerm));


	// calculating principal directions in combinations of the u, v basis

	//umbilic case: kappa1 = kappa2 is taken care of by the a12 = 0

	Vector2D prinCoor1, prinCoor2 ;

	if( a12 == 0 ) {
		if( a11 < a22 ) {
			//DEBUG
			//cout<<"a11 < a22"<<endl;
			//TOCHECK
			prinCoor1 = Vector2D( 0, 1 ) ;
			prinCoor2 = Vector2D( 1, 0 ) ;
		}
		else if( a11 > a22 ) {
			//DEBUG
			//cout<<" a11 > a22"<<endl;
			//TOCHECK
			prinCoor1 = Vector2D( 1, 0 ) ;
			prinCoor2 = Vector2D( 0, 1 ) ;
		}
		else {
			//TOCHECK: seems kappa1 == 0 for this case
			prinCoor1 = Vector2D( 0, 0 ) ;
			prinCoor2 = Vector2D( 0, 0 ) ;
			//DEBUG
			//cout<<"getGradDistance(x+dpx,y,z)"<<getGradDistance(x+dpx,y,z)-getGradDistance(x-dpx,y,z)<<endl;

		}

	}
	else {
		prinCoor1 = Vector2D((a11-a22+sqrt(squaredTerm))/(2*a12),1);

		prinCoor2 = Vector2D(-(-a11+a22+sqrt(squaredTerm))/(2*a12),1);
	}

	// calculating the principal directions

	Vector3D p1 = prinCoor1.getX() * u + prinCoor1.getY() * v;

	Vector3D p2 = prinCoor2.getX() * u + prinCoor2.getY() * v;

	p1.normalize() ;
	p2.normalize() ;

	//DEBUG
	//cout<<"kappa1"<<kappa1<<endl;
	//cout<<"kappa2"<<kappa2<<endl;

	return( Curvature( p1, p2, kappa1, kappa2 ) ) ;

}

bool ImageDistanceMap::getGradCurvature(double x, double y, double z, Vector3D& gradKappa1, Vector3D& gradKappa2) {

	// For calculating the unit distance 

	const Vector3D imageToModelScale(
		1.0/getModelToImageScale()[0],
		1.0/getModelToImageScale()[1],
		1.0/getModelToImageScale()[2] );

	// dp is delta voxel in unit co-ordinates.

	Vector3D dp = imageToModelScale;
	dp = dp / 2.;	// halve the step size as an initial guess
	double dpx = dp.getX();		// faster forms of dp
	double dpy = dp.getY();
	double dpz = dp.getZ();

	// calculating the gradient of gradient of kappa1

	double gradKappa1X =  0.9679 * ( getCurvature(x+dpx, y, z).kappa1 - getCurvature(x-dpx, y, z).kappa1 ) 
		+ 0.0532 * ( getCurvature(x+2*dpx, y, z).kappa1 - getCurvature(x-2*dpx, y, z).kappa1 ) ;
	if(gradKappa1X==0){
         gradKappa1X =  0.9679 * ( getCurvature(x+5*dpx, y, z).kappa1 - getCurvature(x-5*dpx, y, z).kappa1 ) 
		+ 0.0532 * ( getCurvature(x+10*dpx, y, z).kappa1 - getCurvature(x-10*dpx, y, z).kappa1 ) ;
	}


	double gradKappa1Y =  0.9679 * ( getCurvature(x, y+dpy, z).kappa1 - getCurvature(x, y-dpy, z).kappa1 ) 
		+ 0.0532 * ( getCurvature(x, y+2*dpy, z).kappa1 - getCurvature(x, y-2*dpy, z).kappa1 ) ;
	if(gradKappa1Y==0){
	gradKappa1Y =  0.9679 * ( getCurvature(x, y+5*dpy, z).kappa1 - getCurvature(x, y-5*dpy, z).kappa1 ) 
		+ 0.0532 * ( getCurvature(x, y+10*dpy, z).kappa1 - getCurvature(x, y-10*dpy, z).kappa1 ) ;
	}

	double gradKappa1Z =  0.9679 * ( getCurvature(x, y, z+dpz).kappa1 - getCurvature(x, y, z-dpz).kappa1 ) 
		+ 0.0532 * ( getCurvature(x, y, z+2*dpz).kappa1 - getCurvature(x, y, z-2*dpz).kappa1 ) ;
	if(gradKappa1Z==0){
		gradKappa1Z =  0.9679 * ( getCurvature(x, y, z+5*dpz).kappa1 - getCurvature(x, y, z-5*dpz).kappa1 ) 
		+ 0.0532 * ( getCurvature(x, y, z+10*dpz).kappa1 - getCurvature(x, y, z-10*dpz).kappa1 ) ;
	}

	gradKappa1 = Vector3D( gradKappa1X, gradKappa1Y, gradKappa1Z ) ;

	// calculating the gradient of gradient of kappa2

	double gradKappa2X =  0.9679 * ( getCurvature(x+dpx, y, z).kappa2 - getCurvature(x-dpx, y, z).kappa2 ) 
		+ 0.0532 * ( getCurvature(x+2*dpx, y, z).kappa2 - getCurvature(x-2*dpx, y, z).kappa2 ) ;
	if(gradKappa2X==0){
         gradKappa2X =  0.9679 * ( getCurvature(x+5*dpx, y, z).kappa2 - getCurvature(x-5*dpx, y, z).kappa2 ) 
		+ 0.0532 * ( getCurvature(x+10*dpx, y, z).kappa2 - getCurvature(x-10*dpx, y, z).kappa2 ) ;
	}

	double gradKappa2Y =  0.9679 * ( getCurvature(x, y+dpy, z).kappa2 - getCurvature(x, y-dpy, z).kappa2 ) 
		+ 0.0532 * ( getCurvature(x, y+2*dpy, z).kappa2 - getCurvature(x, y-2*dpy, z).kappa2 ) ;
	if(gradKappa2Y==0){
	gradKappa2Y =  0.9679 * ( getCurvature(x, y+5*dpy, z).kappa2 - getCurvature(x, y-5*dpy, z).kappa2 ) 
		+ 0.0532 * ( getCurvature(x, y+10*dpy, z).kappa2 - getCurvature(x, y-10*dpy, z).kappa2 ) ;
	}

	double gradKappa2Z =  0.9679 * ( getCurvature(x, y, z+dpz).kappa2 - getCurvature(x, y, z-dpz).kappa2 ) 
		+ 0.0532 * ( getCurvature(x, y, z+2*dpz).kappa2 - getCurvature(x, y, z-2*dpz).kappa2 ) ;
	if(gradKappa2Z==0){
		gradKappa2Z =  0.9679 * ( getCurvature(x, y, z+5*dpz).kappa2 - getCurvature(x, y, z-5*dpz).kappa2 ) 
		+ 0.0532 * ( getCurvature(x, y, z+10*dpz).kappa2 - getCurvature(x, y, z-10*dpz).kappa2 ) ;
	}

	gradKappa2 = Vector3D( gradKappa2X, gradKappa2Y, gradKappa2Z ) ;

	//DEBUG
	//cout<<"gradKappa1"<<gradKappa1<<endl;
	//cout<<"gradKappa2"<<gradKappa2<<endl;


	return( 1 ) ;

}
