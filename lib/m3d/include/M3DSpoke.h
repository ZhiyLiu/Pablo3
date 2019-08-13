#ifndef _M3D_SPOKE
#define _M3D_SPOKE

#include "Vector3D.h"

class M3DSpoke {

	protected:

		Vector3D X ;	// the position of the tail of the spoke 

		Vector3D U ;	// the unit vector, that is direction of the spoke

		double  r ;		// the length of the spoke

	public:

		M3DSpoke() {
		}

		M3DSpoke( Vector3D _X, Vector3D _U, double _r ) ;

		M3DSpoke( double x, double y, double z, double uX, double uY, double uZ, double _r ) ;

		virtual Vector3D getX() { return X; }

		virtual Vector3D getU() { return U; } 

		virtual double getR() { return r; } 

		virtual Vector3D getS() { return (r * U); } // Returns the scaled spoke direction
 
		virtual Vector3D getB() { return ( X + (r * U) ); } // Returns implied boundary point

		virtual void draw();
        void draw(double r, double g, double b);

        void setX(Vector3D x) {X= x;}
        void setR(double _r) {r = _r;}
        void setU(Vector3D _U) {U = _U;}

        // Zhiyuan Liu: This part is for linking structure of 2D s-rep
        // properties needed in computing linking structure
        int figureID; // the figure this spoke belong to
        int linkTo = -1;   // the figure id this spoke link to, if -1 means unlinked
        int primitiveIndex; // the primitive this spoke attached to
        bool isValid = true;   // some spokes need  to be filtered out before stats analysis
        double linkLength = 1000.0; // linking length
        Vector3D linkPoint; // linking point on external axis
        Vector3D linkPointOnSkeleton; // linking point on another skeleton
        Vector3D linkPointOnBdry;     // linking point on boundary of anther figure
        Vector3D linkPointOnAxis;     // linking point z' extending from another figure, locating on linking axis
        int linkSpokeType; // 0: link to up spoke; 1: link to down spoke, this is spoke type of link
        int spokeType;     // 0: up spoke; 1: down spoke; 2: crest spoke this is spoke type of itself
        int spokeId; // current position in up/down spokes collection
        int linkId; // linked spoke id 
        double linkMinDelta = 10000.0; // the global minimum delta mentioned in Damon's paper.
        void copyFrom(M3DSpoke& s);
};


#endif // _M3D_SPOKE
