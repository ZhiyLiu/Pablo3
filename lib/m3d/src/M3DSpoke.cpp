//#include <GL/glut.h>
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include "M3DSpoke.h"

M3DSpoke::M3DSpoke( Vector3D _X, Vector3D _U, double _r ) {

	X = _X ;
	U = _U ;
	r = _r ;			
}

M3DSpoke::M3DSpoke( double x, double y, double z, double uX, double uY, double uZ, double _r ) {

	X = Vector3D( x, y, z ) ;
	U = Vector3D( uX, uY, uZ ) ;
	r = _r ;

}

void M3DSpoke::draw()
{
	glPushAttrib(GL_CURRENT_BIT);

	glEnable(GL_LIGHTING);
	glPushMatrix();
	glTranslated(X.getX(), X.getY(), X.getZ());

	GLUquadric *qobj = gluNewQuadric();
	glColor3d(1.0,1.0,1.0);
	gluSphere(qobj, r/12.0, 8, 8);
	gluDeleteQuadric(qobj);
	glPopMatrix();
	glDisable(GL_LIGHTING);

	glColor3d(1.0,1.0,1.0);
	glBegin(GL_LINES);
	glVertex3f(X.getX(), X.getY(), X.getZ());
	glVertex3f(this->getB().getX(), this->getB().getY(), this->getB().getZ());
	glEnd();

}

void M3DSpoke::draw(double r, double g, double b)
{
    
	// glPushAttrib(GL_CURRENT_BIT);

	// glEnable(GL_LIGHTING);
	// glPushMatrix();
	// glTranslated(X.getX(), X.getY(), X.getZ());

	// GLUquadric *qobj = gluNewQuadric();
	// glColor3d(1.0,1.0,1.0);
	// gluSphere(qobj, r/12.0, 8, 8);
	// gluDeleteQuadric(qobj);
	// glPopMatrix();
	// glDisable(GL_LIGHTING);

	glColor3d(r, g, b);
	glBegin(GL_LINES);
	glVertex3f(X.getX(), X.getY(), X.getZ());
	glVertex3f(this->getB().getX(), this->getB().getY(), this->getB().getZ());
	glEnd();
}

void M3DSpoke::copyFrom(M3DSpoke& s)
{
    setX(s.getX());
    setR(s.getR());
    setU(s.getU());
    linkTo = s.linkTo;
    figureID = s.figureID;
    primitiveIndex = s.primitiveIndex;
    linkLength = s.linkLength;
    linkPoint = s.linkPoint;
}

