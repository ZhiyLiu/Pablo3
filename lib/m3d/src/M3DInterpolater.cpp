#include <math.h>
#include "M3DInterpolater.h"
//#include <GL/gl.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_4.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_4x4.h>
#include <vnl/vnl_double_2x2.h>
#include <vnl/vnl_double_2x3.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_1x3.h>
#include <vnl/vnl_double_3x1.h>
#include <vnl/vnl_double_1x1.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vtkLine.h>
#include <stdio.h>

M3DInterpolater::M3DInterpolater()
{
	_figure = NULL;
}

M3DInterpolater::M3DInterpolater ( M3DFigure *figure )
{
    _figure = figure;
}

void M3DInterpolater::convertVNLToVector3D ( vnl_double_3x1 input, Vector3D &out )
{
    out.setX(input(0,0));
    out.setY(input(1,0));
    out.setZ(input(2,0));
}
void M3DInterpolater::convertVector3DToVNL ( Vector3D vec, vnl_double_3x1 &out )
{
    out(0,0) = vec.getX();
    out(1,0) = vec.getY();
    out(2,0) = vec.getZ();
}

vnl_double_3x1 M3DInterpolater::convertVNL4to3(vnl_matrix_fixed<double,4,1> q)
{
	vnl_double_3x1 result;

	if ( fabs(q(0,0)) > 1e-10 )
		std::cout << "Warning: possible illegal conversion from 4 to 3" << std::endl;

	result(0,0) = q(1,0);  result(1,0) = q(2,0);  result(2,0) = q(3,0);

	return result;
}
vnl_matrix_fixed<double,4,1> M3DInterpolater::convertVNL3to4(vnl_double_3x1 q)
{
	vnl_matrix_fixed<double,4,1> result;

	result(0,0) = 0;  result(1,0) = q(0,0);  result(2,0) = q(1,0);  result(3,0) = q(2,0);

	return result;
}

void M3DInterpolater::displayInterpolatedBoundary ( M3DFigure *figure, int level )
{
//    glBegin ( GL_QUADS );
//    glNormal3d ( 0.0,0.0,1.0 );
//    glVertex3d ( 0.0,0.0,0.0 );
//    glVertex3d ( 1.0,0.0,0.0 );
//    glVertex3d ( 0.0,1.0,0.0 );
//    glVertex3d ( 1.0,1.0,0.0 );
//    glEnd();
}

void M3DInterpolater::qlog(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> &result)
{
	double theta = acos(q(0,0));

	result(0,0) = 0;
	result(1,0) = (theta*q(1,0))/sin(theta);
	result(2,0) = (theta*q(2,0))/sin(theta);
	result(3,0) = (theta*q(3,0))/sin(theta);
}

void M3DInterpolater::qexp(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> &result)
{
	double a = exp(q(0,0));
	vnl_double_3x1 v;
	v(0,0) = q(1,0);
	v(1,0) = q(2,0);
	v(2,0) = q(3,0);
	double nv = v.get_column(0).magnitude();

	result(0,0) = a * cos(nv);
	result(1,0) = a * ( sin(nv) * ( q(1,0)/nv ) );
	result(2,0) = a * ( sin(nv) * ( q(2,0)/nv ) );
	result(3,0) = a * ( sin(nv) * ( q(3,0)/nv ) );
}

void M3DInterpolater::qpow(vnl_matrix_fixed<double,4,1> q, double t, vnl_matrix_fixed<double,4,1> &result)
{
	vnl_matrix_fixed<double,4,1> result2;
	qlog(q,result2);
	qexp(t*result2,result);
}

void M3DInterpolater::qmul(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> p, vnl_matrix_fixed<double,4,1> &result)
{
	result(0,0) = q(0,0)*p(0,0) - q(1,0)*p(1,0) - q(2,0)*p(2,0) - q(3,0)*p(3,0);
	result(1,0) = q(0,0)*p(1,0) + q(1,0)*p(0,0) + q(2,0)*p(3,0) - q(3,0)*p(2,0);
	result(2,0) = q(0,0)*p(2,0) - q(1,0)*p(3,0) + q(2,0)*p(0,0) + q(3,0)*p(1,0);
	result(3,0) = q(0,0)*p(3,0) + q(1,0)*p(2,0) - q(2,0)*p(1,0) + q(3,0)*p(0,0);
}

void M3DInterpolater::qcon(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> &result)
{
	result(0,0) =  q(0,0);
	result(1,0) = -q(1,0);
	result(2,0) = -q(2,0);
	result(3,0) = -q(3,0);
}

vnl_matrix_fixed<double,4,1> M3DInterpolater::qlog(vnl_matrix_fixed<double,4,1> q)
{
	vnl_matrix_fixed<double,4,1> result;

	double theta = acos(q(0,0));

	result(0,0) = 0;
	result(1,0) = (theta*q(1,0))/sin(theta);
	result(2,0) = (theta*q(2,0))/sin(theta);
	result(3,0) = (theta*q(3,0))/sin(theta);

	return result;
}
vnl_matrix_fixed<double,4,1> M3DInterpolater::qpow(vnl_matrix_fixed<double,4,1> q, double t)
{
	vnl_matrix_fixed<double,4,1> result;

	result = qexp(t*qlog(q));

	return result;
}
vnl_matrix_fixed<double,4,1> M3DInterpolater::qexp(vnl_matrix_fixed<double,4,1> q)
{
	vnl_matrix_fixed<double,4,1> result;
	double q0 = q(0,0); double q1 = q(1,0); double q2 = q(2,0); double q3 = q(3,0);

	double a = exp(q0);
	double nv = sqrt(q1*q1 + q2*q2 + q3*q3);
	double snv = sin(nv);

	result(0,0) = a * cos(nv);
	result(1,0) = a * ( snv * ( q1/nv ) );
	result(2,0) = a * ( snv * ( q2/nv ) );
	result(3,0) = a * ( snv * ( q3/nv ) );

	return result;
}
vnl_matrix_fixed<double,4,1> M3DInterpolater::qmul(vnl_matrix_fixed<double,4,1> q, vnl_matrix_fixed<double,4,1> p)
{
	vnl_matrix_fixed<double,4,1> result;
	double q0 = q(0,0); double q1 = q(1,0); double q2 = q(2,0); double q3 = q(3,0);
	double p0 = p(0,0); double p1 = p(1,0); double p2 = p(2,0); double p3 = p(3,0);

	result(0,0) = q0*p0 - q1*p1 - q2*p2 - q3*p3;
	result(1,0) = q0*p1 + q1*p0 + q2*p3 - q3*p2;
	result(2,0) = q0*p2 - q1*p3 + q2*p0 + q3*p1;
	result(3,0) = q0*p3 + q1*p2 - q2*p1 + q3*p0;

	return result;
}
vnl_matrix_fixed<double,4,1> M3DInterpolater::qcon(vnl_matrix_fixed<double,4,1> q)
{
	vnl_matrix_fixed<double,4,1> result;

	result(0,0) =  q(0,0);
	result(1,0) = -q(1,0);
	result(2,0) = -q(2,0);
	result(3,0) = -q(3,0);

	return result;
}

vnl_matrix_fixed<double,4,1> M3DInterpolater::qslerp(vnl_matrix_fixed<double,4,1> p, vnl_matrix_fixed<double,4,1> q, double h)
{
	vnl_matrix_fixed<double,4,1> result;

	result = qmul(p,qpow(qmul(qcon(p),q),h));

	return result;
}

void M3DInterpolater::slerp(vnl_double_3x1 U1, vnl_double_3x1 U2, double u, vnl_double_3x1 &result)
{
	double phi = acos((U1.transpose()*U2)(0,0));

	result = ( sin((1-u)*phi)/sin(phi) )*U1 + ( sin(u*phi)/sin(phi) )*U2;
}

void M3DInterpolater::bislerp(vnl_double_3x1 U11, vnl_double_3x1 U12, vnl_double_3x1 U21, vnl_double_3x1 U22, double u, double v, vnl_double_3x1 &result)
{
	// First, do two slerps in u direction, then one between those in v direction

	vnl_double_3x1 tempU1, tempU2;

	slerp(U11, U21, u, tempU1);
	slerp(U12, U22, u, tempU2);
	slerp(tempU1, tempU2, v, result);

}

vnl_double_3x1 M3DInterpolater::squad(vnl_double_3x1 u0, vnl_double_3x1 u1, vnl_double_3x1 u2, vnl_double_3x1 u3, double h)
{
	vnl_matrix_fixed<double,4,1> q0, q1, q2, q3;
	q0 = convertVNL3to4(u0);
	q1 = convertVNL3to4(u1);
	q2 = convertVNL3to4(u2);
	q3 = convertVNL3to4(u3);

	vnl_matrix_fixed<double,4,1> si, sip;
	si = qmul( q1, qexp( -0.25* ( qlog(qmul(qcon(q1),q2)) + qlog(qmul(qcon(q1),q0)) ) ) );
	sip = qmul( q2, qexp( -0.25* ( qlog(qmul(qcon(q2),q3)) + qlog(qmul(qcon(q2),q1)) ) ) );

	vnl_double_3x1 ssi, ssip;
	ssi = convertVNL4to3(si);
	ssip = convertVNL4to3(sip);

	vnl_double_3x1 test1, test2, result;
	slerp(u1, u2, h, test1);
	slerp(ssi, ssip, h, test2);
	slerp(test1, test2, 2*h*(1-h), result);

	return result;
}

vnl_double_3x1 M3DInterpolater::dsquad(vnl_double_3x1 u0, vnl_double_3x1 u1, vnl_double_3x1 u2, vnl_double_3x1 u3, double h)
{
	vnl_matrix_fixed<double,4,1> q0, q1, q2, q3;
	q0 = convertVNL3to4(u0);
	q1 = convertVNL3to4(u1);
	q2 = convertVNL3to4(u2);
	q3 = convertVNL3to4(u3);

	vnl_matrix_fixed<double,4,1> si, sip;
	si = qmul( q1, qexp( -0.25* ( qlog(qmul(qcon(q1),q2)) + qlog(qmul(qcon(q1),q0)) ) ) );
	sip = qmul( q2, qexp( -0.25* ( qlog(qmul(qcon(q2),q3)) + qlog(qmul(qcon(q2),q1)) ) ) );

	vnl_matrix_fixed<double,4,1> term1, term2, term3, term4, term5;
	vnl_double_3x1 temp1,temp2;

	double hp = h + 1e-5;
	double hm = h - 1e-5;

	slerp(u1,u2,h,temp1);
	term1 = convertVNL3to4(temp1);

	term2 = qlog( qmul( qcon(q1), q2 ) );

	term3 = qpow( g(q1,q2,si,sip,h), 2*h*(1-h) );

	term4 = convertVNL3to4(temp1);

	term5 = 1e5 * ( qpow( g(q1,q2,si,sip,hp), 2*hp*(1-hp) ) - qpow( g(q1,q2,si,sip,hm), 2*hm*(1-hm) ) );

	vnl_matrix_fixed<double,4,1> res = qmul( qmul(term1,term2), term3 ) + qmul(term4,term5);

	vnl_double_3x1 result = convertVNL4to3(res);

	return result;

}

vnl_matrix_fixed<double,4,1> M3DInterpolater::g(vnl_matrix_fixed<double,4,1> qi, vnl_matrix_fixed<double,4,1> qip,
		vnl_matrix_fixed<double,4,1> si, vnl_matrix_fixed<double,4,1> sip, double h)
{
	vnl_double_3x1 qqi, qqip, ssi, ssip;
	qqi = convertVNL4to3(qi);
	qqip = convertVNL4to3(qip);
	ssi = convertVNL4to3(si);
	ssip = convertVNL4to3(sip);

	vnl_double_3x1 temp1, temp2;
	slerp(qqi,qqip,h,temp1);
	slerp(ssi,ssip,h,temp2);

	vnl_matrix_fixed<double,4,1> result, temp3, temp4;

	temp3 = convertVNL3to4(temp1);
	temp4 = convertVNL3to4(temp2);

	result = qmul( qcon(temp3), temp4 );

	return result;
}

vnl_matrix_fixed<double,4,1> M3DInterpolater::dg(vnl_matrix_fixed<double,4,1> qi, vnl_matrix_fixed<double,4,1> qip,
		vnl_matrix_fixed<double,4,1> si, vnl_matrix_fixed<double,4,1> sip, double h)
{
	vnl_matrix_fixed<double,4,1> gp = g(qi,qip,si,sip,h+1e-5);
	vnl_matrix_fixed<double,4,1> gm = g(qi,qip,si,sip,h-1e-5);

	return 1e5*(gp-gm);
}


double  M3DInterpolater::getR(M3DPrimitive* p, int side)
{
    if(side == 0)
    {
        return dynamic_cast<M3DQuadPrimitive*>(p)->getR0();
    }
    else
    {
         return dynamic_cast<M3DQuadPrimitive*>(p)->getR1();
    }

}

void M3DInterpolater::getU(M3DPrimitive* p, int side, vnl_double_3x1 &u)
{
    if(side ==0 )
    {
        convertVector3DToVNL(dynamic_cast<M3DQuadPrimitive*>(p)->getU0(), u);
    }
    else
    {
        convertVector3DToVNL(dynamic_cast<M3DQuadPrimitive*>(p)->getU1(), u);
    }

}
void M3DInterpolater::setR(M3DPrimitive* p, int side, double r)
{
    if(side == 0)
    {
        dynamic_cast<M3DQuadPrimitive*>(p)->setR0(r);
    }
    else
    {
         dynamic_cast<M3DQuadPrimitive*>(p)->setR1(r);
    }

}

void M3DInterpolater::setU(M3DPrimitive* p, int side, vnl_double_3x1 u)
{
    Vector3D vecU(u(0,0), u(1,0), u(2,0));
    if(side == 0)
    {
        dynamic_cast<M3DQuadPrimitive*>(p)->setU0(vecU);
    }
    else
    {
        dynamic_cast<M3DQuadPrimitive*>(p)->setU1(vecU);
    }

}
void M3DInterpolater::interpolateX(M3DFigure *figure,  int uBase, int vBase, double u, double v, int side, vnl_double_3x1& output)
{
    M3DQuadFigure *fig = dynamic_cast<M3DQuadFigure*>(figure);
    vnl_double_3x1 x11, x12, x21, x22;

	convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase)->getX(), x11);
	convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase+1)->getX(), x12);
	convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase)->getX(), x21);
	convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase+1)->getX(), x22);

	vnl_double_3x1 dxdu11, dxdu12, dxdu21, dxdu22;
	vnl_double_3x1 dxdv11, dxdv12, dxdv21, dxdv22;

	convertVector3DToVNL(fig->computeDXdu(uBase,vBase), dxdu11);
	convertVector3DToVNL(fig->computeDXdu(uBase,vBase+1), dxdu12);
	convertVector3DToVNL(fig->computeDXdu(uBase+1,vBase), dxdu21);
	convertVector3DToVNL(fig->computeDXdu(uBase+1,vBase+1), dxdu22);

	convertVector3DToVNL(fig->computeDXdv(uBase,vBase), dxdv11);
	convertVector3DToVNL(fig->computeDXdv(uBase,vBase+1), dxdv12);
	convertVector3DToVNL(fig->computeDXdv(uBase+1,vBase), dxdv21);
	convertVector3DToVNL(fig->computeDXdv(uBase+1,vBase+1), dxdv22);

    vnl_double_4x4 hx, hy, hz;

	hx(0,0) = x11(0,0); 	hx(0,1) = x12(0,0);
	hx(1,0) = x21(0,0);	 	hx(1,1) = x22(0,0);
	hx(2,0) = dxdu11(0,0); 	hx(2,1) = dxdu12(0,0);
	hx(3,0) = dxdu21(0,0); 	hx(3,1) = dxdu22(0,0);
	hx(0,2) = dxdv11(0,0);	hx(0,3) = dxdv12(0,0);
	hx(1,2) = dxdv21(0,0);	hx(1,3) = dxdv22(0,0);
	hx(2,2) = 0; hx(2,3) = 0; hx(3,2) = 0; hx(3,3) = 0;


	hy(0,0) = x11(1,0);		hy(0,1) = x12(1,0);
	hy(1,0) = x21(1,0);	 	hy(1,1) = x22(1,0);
	hy(2,0) = dxdu11(1,0); 	hy(2,1) = dxdu12(1,0);
	hy(3,0) = dxdu21(1,0); 	hy(3,1) = dxdu22(1,0);
	hy(0,2) = dxdv11(1,0);	hy(0,3) = dxdv12(1,0);
	hy(1,2) = dxdv21(1,0);	hy(1,3) = dxdv22(1,0);
	hy(2,2) = 0; hy(2,3) = 0; hy(3,2) = 0; hy(3,3) = 0;

	hz(0,0) = x11(2,0);		hz(0,1) = x12(2,0);
	hz(1,0) = x21(2,0);	 	hz(1,1) = x22(2,0);
	hz(2,0) = dxdu11(2,0); 	hz(2,1) = dxdu12(2,0);
	hz(3,0) = dxdu21(2,0); 	hz(3,1) = dxdu22(2,0);
	hz(0,2) = dxdv11(2,0);	hz(0,3) = dxdv12(2,0);
	hz(1,2) = dxdv21(2,0);	hz(1,3) = dxdv22(2,0);
	hz(2,2) = 0; hz(2,3) = 0; hz(3,2) = 0; hz(3,3) = 0;

	vnl_matrix_fixed<double,4,1> hu, hv;
	hu(0,0) = h1(u);  hu(1,0) = h2(u);  hu(2,0) = h3(u);  hu(3,0) = h4(u);
	hv(0,0) = h1(v);  hv(1,0) = h2(v);  hv(2,0) = h3(v);  hv(3,0) = h4(v);

	vnl_double_1x1 xn = hu.transpose() * hx * hv;
	vnl_double_1x1 yn = hu.transpose() * hy * hv;
	vnl_double_1x1 zn = hu.transpose() * hz * hv;

    output(0,0) = xn(0,0);
    output(1,0) = yn(0,0);
    output(2,0) = zn(0,0);
}

void M3DInterpolater::compute2ndDerivative(vnl_double_3x1 startU, vnl_double_3x1 endU, vnl_double_3x1 targetU, double d, vnl_double_3x1& out)
{
    double del = 1e-5;
    vnl_double_3x1 Upv1, Upv2, Upv3, Upv4, Upv5;
    slerp(startU, endU, d+2*del, Upv1);
    slerp(startU, endU, d+del, Upv2);
    //slerp(startU, endU, d, Upv3);
    slerp(startU, endU, d-del, Upv4);
    slerp(startU, endU, d-2*del, Upv5);

    out = 0.25 * (Upv5 + Upv1 - 2.0 * targetU);

}
// d is distance from end to start
double M3DInterpolater::computeMiddleS(vnl_double_3x1 startU, vnl_double_3x1 endU, vnl_double_3x1 startS, vnl_double_3x1 endS, double d, vnl_double_3x1& uMiddle)
{

    // 1. at first d is the length of interval
    vnl_double_3x1 Uvv_start, Uvv_end;
    compute2ndDerivative(startU, endU, endU, d, Uvv_end);
    compute2ndDerivative(startU, endU, startU, 0, Uvv_start);

    // 2. then d becomes distance from middle point to left-most point
    d /= 2;
    slerp(startU, endU, d, uMiddle);
    vnl_double_3x1 avg = 0.5 * (startS + endS);
    vnl_double_1x1 result = uMiddle.transpose() * avg - d*d*0.25*(startS.transpose() * Uvv_start + endS.transpose() * Uvv_end);

    return result(0,0);
}

// recursive divide into 4 share
// lamda is power function of 1/2, first invocation is 1, equally space between u dir and v dir
double M3DInterpolater::interpolateS(M3DFigure *figure, M3DPrimitive* Sp11, M3DPrimitive* Sp12, M3DPrimitive* Sp21, M3DPrimitive* Sp22, double u, double v, int side, double lamda, Vector3D& interpolatedU)
{
	M3DQuadFigure *fig = dynamic_cast<M3DQuadFigure*>(figure);

    // compute r(p+dv)
    vnl_double_3x1 U11, U12, U21, U22;
    vnl_double_3x1 avgTop, avgLeft, avgRight, avgBottom;
    double r11, r12, r21, r22;

    r11 = getR(Sp11, side);
    r12 = getR(Sp12, side);
    r21 = getR(Sp21, side);
    r22 = getR(Sp22, side);

    getU(Sp11, side, U11);
    getU(Sp12, side, U12);
    getU(Sp21, side, U21);
    getU(Sp22, side, U22);

    // middle spoke on top edge
    vnl_double_3x1 uMidTop;
    double rMidTop = computeMiddleS(U11, U12, r11 * U11, r12 * U12, lamda, uMidTop);

    // middle spoke on left edge
    vnl_double_3x1 uMidLeft;
    double rMidLeft = computeMiddleS(U11, U21, r11*U11, r21 * U21, lamda, uMidLeft);

    // middle spoke on right edge
    vnl_double_3x1 uMidRight;
    double rMidRight = computeMiddleS(U12, U22, r12 * U12, r22 * U22, lamda, uMidRight);

    // middle spoke on bottom edge
    vnl_double_3x1 uMidBot;
    double rMidBot = computeMiddleS(U21, U22, r21 * U21, r22 * U22, lamda, uMidBot);

    //center spoke in this quad
    vnl_double_3x1 uCenterA, uCenterB;
    double rCenterA = computeMiddleS(uMidTop, uMidBot, rMidTop*uMidTop, rMidBot*uMidBot, lamda, uCenterA);
    double rCenterB = computeMiddleS(uMidLeft, uMidRight, rMidLeft*uMidLeft, rMidRight*uMidRight, lamda, uCenterB);
    double rCenter = 0.5 * (rCenterA + rCenterB);
    vnl_double_3x1 uCenter = 0.5 * (uCenterA + uCenterB);

    double halfDist = lamda / 2;
    double tolerance = 1e-6;
    if(abs(u - halfDist ) <= tolerance && abs(v - halfDist) <= tolerance)
    {
        // if the target spoke locates in the center of this quad, return center spoke
        convertVNLToVector3D(uCenter, interpolatedU);
        return rCenter;
    }
    else if(abs(u) < tolerance && abs(v) < tolerance)
    {
        // left top corner
        convertVNLToVector3D(U11, interpolatedU);
        return r11;
    }
    else if(abs(u-lamda) < tolerance && abs(v) < tolerance)
    {
        // left bottom corner
        convertVNLToVector3D(U21, interpolatedU);
        return r21;
    }
    else if(abs(v-lamda) < tolerance && abs(u) < tolerance)
    {
        // right top corner
        convertVNLToVector3D(U12, interpolatedU);
        return r12;
    }
    else if(abs(v-lamda) < tolerance && abs(u-lamda) < tolerance)
    {
        // right bottom corner
        convertVNLToVector3D(U22, interpolatedU);
        return r22;
    }
    else if(abs(u - halfDist) <= tolerance && abs(v) < tolerance)
    {
        // on the left edge
        convertVNLToVector3D(uMidLeft, interpolatedU);
        return rMidLeft;
    }
    else if(abs(u) < tolerance && abs(v - halfDist) < tolerance)
    {
        // on the top edge
        convertVNLToVector3D(uMidTop, interpolatedU);
        return rMidTop;
    }
    else if(abs(u - halfDist) < tolerance && abs(v-lamda) < tolerance)
    {
        // right edge
        convertVNLToVector3D(uMidRight, interpolatedU);
        return rMidRight;
    }
    else if(abs(u-lamda) < tolerance && abs(v-halfDist) < tolerance)
    {
        // bot edge
        convertVNLToVector3D(uMidBot, interpolatedU);
        return rMidBot;
    }
    else
    {
        // construct primitives
        M3DQuadPrimitive midTop, midLeft, midRight, midBot, center;
        setR(&midTop, side, rMidTop);
        setR(&midLeft, side, rMidLeft);
        setR(&midRight, side, rMidRight);
        setR(&midBot, side, rMidBot);
        setR(&center, side, rCenter);

        setU(&midTop, side, uMidTop);
        setU(&midLeft, side, uMidLeft);
        setU(&midRight, side, uMidRight);
        setU(&midBot, side, uMidBot);
        setU(&center, side, uCenter);

        // if the target spoke locates in the 1st quadrant
        if(u < halfDist && v > halfDist)
        {
            return interpolateS(figure, &midTop, Sp12, &center, &midRight, u, v-halfDist, side, lamda /2, interpolatedU);
        }

        // if the target spoke locates in the 2nd quadrant
        else if(u < halfDist && v < halfDist)
        {
            return interpolateS(figure, Sp11, &midTop, &midLeft, &center, u, v, side, lamda / 2, interpolatedU);

        }

        // if the target spoke locates in the 3rd quadrant
        else if(u > halfDist && v < halfDist)
        {
            return interpolateS(figure, &midLeft, &center, Sp21, &midBot, u - halfDist, v, side, lamda/2, interpolatedU);

        }
        // if the target spoke locates in the 4th quadrant
        else if(u > halfDist && v > halfDist)
        {
            return interpolateS(figure, &center, &midRight, &midBot, Sp22, u - halfDist, v - halfDist, side, lamda/2, interpolatedU);
        }
    }

}

// This is implementation of math derivated by Steve and Zhiyuan in 2018
M3DSpoke M3DInterpolater::interpolate(M3DFigure *figure, double u, double v, int side)
{
    // interpolate r(u,v)
    std::cout << "U: " << u << " " << v << " " << side << std::endl;
	int uBase = (int) floor(u);
	int vBase = (int) floor(v);

	u = u - uBase;
	v = v - vBase;

	M3DQuadFigure *fig = dynamic_cast<M3DQuadFigure*>(figure);

	if (uBase == fig->getRowCount() - 1)
	{
		uBase = uBase - 1;
		u = 1;
	}

	if (vBase == fig->getColumnCount() - 1)
	{
		vBase = vBase - 1;
		v = 1;
	}

    // compute r(p+dv)
    M3DPrimitive *Sp11 = fig->getPrimitivePtr(uBase,vBase);
    M3DPrimitive *Sp12 = fig->getPrimitivePtr(uBase, vBase+1);
    M3DPrimitive *Sp21 = fig->getPrimitivePtr(uBase+1, vBase);
    M3DPrimitive *Sp22 = fig->getPrimitivePtr(uBase+1, vBase+1);

    Vector3D interpolatedU;
    double r_uv = interpolateS(figure, Sp11, Sp12, Sp21, Sp22, u, v, side, 1, interpolatedU);
    
    // Now we have r and U of this spoke, we need interpolate its skeletal point
    vnl_double_3x1 pt;
    interpolateX(fig, uBase, vBase, u, v, side, pt);
    M3DSpoke retSpoke(pt(0,0),pt(1,0),pt(2,0), interpolatedU.getX(), interpolatedU.getY(), interpolatedU.getZ(), r_uv);
    
    return retSpoke;
}

// This interpolation comes from Jared's dissertation. The result is good. Yet the math is wrong.
M3DSpoke M3DInterpolater::interpolateSpoke_new( M3DFigure *figure, double u, double v, int side )
{
	std::cout << "U: " << u << " " << v << " " << side << std::endl;
	int uBase = (int) floor(u);
	int vBase = (int) floor(v);

	u = u - uBase;
	v = v - vBase;

	M3DQuadFigure *fig = dynamic_cast<M3DQuadFigure*>(figure);

	if (uBase == fig->getRowCount() - 1)
	{
		uBase = uBase - 1;
		u = 1;
	}

	if (vBase == fig->getColumnCount() - 1)
	{
		vBase = vBase - 1;
		v = 1;
	}
	
	bool onEdge = (uBase == 0) || (vBase == 0) ||
				  (uBase == fig->getRowCount() - 2) || 
				  (vBase == fig->getColumnCount() - 2);

	vnl_double_3x1 x11, x12, x21, x22;

	convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase)->getX(), x11);
	convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase+1)->getX(), x12);
	convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase)->getX(), x21);
	convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase+1)->getX(), x22);

	vnl_double_3x1 dxdu11, dxdu12, dxdu21, dxdu22;
	vnl_double_3x1 dxdv11, dxdv12, dxdv21, dxdv22;

	convertVector3DToVNL(fig->computeDXdu(uBase,vBase), dxdu11);
	convertVector3DToVNL(fig->computeDXdu(uBase,vBase+1), dxdu12);
	convertVector3DToVNL(fig->computeDXdu(uBase+1,vBase), dxdu21);
	convertVector3DToVNL(fig->computeDXdu(uBase+1,vBase+1), dxdu22);

	convertVector3DToVNL(fig->computeDXdv(uBase,vBase), dxdv11);
	convertVector3DToVNL(fig->computeDXdv(uBase,vBase+1), dxdv12);
	convertVector3DToVNL(fig->computeDXdv(uBase+1,vBase), dxdv21);
	convertVector3DToVNL(fig->computeDXdv(uBase+1,vBase+1), dxdv22);

	vnl_double_3x1 U11, U12, U21, U22;
	double r11, r12, r21, r22;
	double drdu11, drdu12, drdu21, drdu22;
	double drdv11, drdv12, drdv21, drdv22;

	vnl_double_3x1 dUdu11, dUdu12, dUdu21, dUdu22;
	vnl_double_3x1 dUdv11, dUdv12, dUdv21, dUdv22;
	
	vnl_double_3x1 U00, U01, U02, U03, U10, U13, U20, U23, U30, U31, U32, U33;

	if (side == 0)
	{
		convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase)->getU0(),U11);
		convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase+1)->getU0(),U12);
		convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase)->getU0(),U21);
		convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase+1)->getU0(),U22);

		r11 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase,vBase))->getR0();
		r12 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase,vBase+1))->getR0();
		r21 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase+1,vBase))->getR0();
		r22 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase+1,vBase+1))->getR0();

		drdu11 = fig->computeDR0du(uBase,vBase);
		drdu12 = fig->computeDR0du(uBase,vBase+1);
		drdu21 = fig->computeDR0du(uBase+1,vBase);
		drdu22 = fig->computeDR0du(uBase+1,vBase+1);

		drdv11 = fig->computeDR0dv(uBase,vBase);
		drdv12 = fig->computeDR0dv(uBase,vBase+1);
		drdv21 = fig->computeDR0dv(uBase+1,vBase);
		drdv22 = fig->computeDR0dv(uBase+1,vBase+1);

		convertVector3DToVNL(fig->computeDU0du(uBase,vBase),dUdu11);
		convertVector3DToVNL(fig->computeDU0du(uBase,vBase+1),dUdu12);
		convertVector3DToVNL(fig->computeDU0du(uBase+1,vBase),dUdu21);
		convertVector3DToVNL(fig->computeDU0du(uBase+1,vBase+1),dUdu22);

		convertVector3DToVNL(fig->computeDU0dv(uBase,vBase),dUdv11);
		convertVector3DToVNL(fig->computeDU0dv(uBase,vBase+1),dUdv12);
		convertVector3DToVNL(fig->computeDU0dv(uBase+1,vBase),dUdv21);
		convertVector3DToVNL(fig->computeDU0dv(uBase+1,vBase+1),dUdv22);
		
		if (!onEdge)
		{
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase-1)->getU0(), U00);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase)->getU0(), U01);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase+1)->getU0(), U02);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase+2)->getU0(), U03);
			
			convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase-1)->getU0(), U10);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase+2)->getU0(), U13);
			
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase-1)->getU0(), U20);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase+2)->getU0(), U23);
			
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase-1)->getU0(), U30);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase)->getU0(), U31);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase+1)->getU0(), U32);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase+2)->getU0(), U33);
		}

	}
	else
	{
		convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase)->getU1(),U11);
		convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase+1)->getU1(),U12);
		convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase)->getU1(),U21);
		convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase+1)->getU1(),U22);

		r11 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase,vBase))->getR1();
		r12 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase,vBase+1))->getR1();
		r21 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase+1,vBase))->getR1();
		r22 = dynamic_cast<M3DQuadPrimitive*>(fig->getPrimitivePtr(uBase+1,vBase+1))->getR1();

		drdu11 = fig->computeDR1du(uBase,vBase);
		drdu12 = fig->computeDR1du(uBase,vBase+1);
		drdu21 = fig->computeDR1du(uBase+1,vBase);
		drdu22 = fig->computeDR1du(uBase+1,vBase+1);

		drdv11 = fig->computeDR1dv(uBase,vBase);
		drdv12 = fig->computeDR1dv(uBase,vBase+1);
		drdv21 = fig->computeDR1dv(uBase+1,vBase);
		drdv22 = fig->computeDR1dv(uBase+1,vBase+1);

		convertVector3DToVNL(fig->computeDU1du(uBase,vBase),dUdu11);
		convertVector3DToVNL(fig->computeDU1du(uBase,vBase+1),dUdu12);
		convertVector3DToVNL(fig->computeDU1du(uBase+1,vBase),dUdu21);
		convertVector3DToVNL(fig->computeDU1du(uBase+1,vBase+1),dUdu22);

		convertVector3DToVNL(fig->computeDU1dv(uBase,vBase),dUdv11);
		convertVector3DToVNL(fig->computeDU1dv(uBase,vBase+1),dUdv12);
		convertVector3DToVNL(fig->computeDU1dv(uBase+1,vBase),dUdv21);
		convertVector3DToVNL(fig->computeDU1dv(uBase+1,vBase+1),dUdv22);
		
		if (!onEdge)
		{
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase-1)->getU1(), U00);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase)->getU1(), U01);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase+1)->getU1(), U02);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase-1,vBase+2)->getU1(), U03);
			
			convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase-1)->getU1(), U10);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase,vBase+2)->getU1(), U13);
			
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase-1)->getU1(), U20);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+1,vBase+2)->getU1(), U23);
			
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase-1)->getU1(), U30);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase)->getU1(), U31);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase+1)->getU1(), U32);
			convertVector3DToVNL(fig->getPrimitivePtr(uBase+2,vBase+2)->getU1(), U33);
		}
	}

	vnl_double_4x4 hx, hy, hz;

	hx(0,0) = x11(0,0); 	hx(0,1) = x12(0,0);
	hx(1,0) = x21(0,0);	 	hx(1,1) = x22(0,0);
	hx(2,0) = dxdu11(0,0); 	hx(2,1) = dxdu12(0,0);
	hx(3,0) = dxdu21(0,0); 	hx(3,1) = dxdu22(0,0);
	hx(0,2) = dxdv11(0,0);	hx(0,3) = dxdv12(0,0);
	hx(1,2) = dxdv21(0,0);	hx(1,3) = dxdv22(0,0);
	hx(2,2) = 0; hx(2,3) = 0; hx(3,2) = 0; hx(3,3) = 0;


	hy(0,0) = x11(1,0);		hy(0,1) = x12(1,0);
	hy(1,0) = x21(1,0);	 	hy(1,1) = x22(1,0);
	hy(2,0) = dxdu11(1,0); 	hy(2,1) = dxdu12(1,0);
	hy(3,0) = dxdu21(1,0); 	hy(3,1) = dxdu22(1,0);
	hy(0,2) = dxdv11(1,0);	hy(0,3) = dxdv12(1,0);
	hy(1,2) = dxdv21(1,0);	hy(1,3) = dxdv22(1,0);
	hy(2,2) = 0; hy(2,3) = 0; hy(3,2) = 0; hy(3,3) = 0;

	hz(0,0) = x11(2,0);		hz(0,1) = x12(2,0);
	hz(1,0) = x21(2,0);	 	hz(1,1) = x22(2,0);
	hz(2,0) = dxdu11(2,0); 	hz(2,1) = dxdu12(2,0);
	hz(3,0) = dxdu21(2,0); 	hz(3,1) = dxdu22(2,0);
	hz(0,2) = dxdv11(2,0);	hz(0,3) = dxdv12(2,0);
	hz(1,2) = dxdv21(2,0);	hz(1,3) = dxdv22(2,0);
	hz(2,2) = 0; hz(2,3) = 0; hz(3,2) = 0; hz(3,3) = 0;

	vnl_matrix_fixed<double,4,1> hu, hv;
	hu(0,0) = h1(u);  hu(1,0) = h2(u);  hu(2,0) = h3(u);  hu(3,0) = h4(u);
	hv(0,0) = h1(v);  hv(1,0) = h2(v);  hv(2,0) = h3(v);  hv(3,0) = h4(v);

	vnl_double_1x1 xn = hu.transpose() * hx * hv;
	vnl_double_1x1 yn = hu.transpose() * hy * hv;
	vnl_double_1x1 zn = hu.transpose() * hz * hv;

	//Vector3D newX = Vector3D(xn(0,0),yn(0,0), zn(0,0));

	vnl_double_4x4 ux, uy, uz;

	ux(0, 0) = U11(0, 0);		ux(0, 1) = U12(0, 0);
	ux(1, 0) = U21(0, 0);	 	ux(1, 1) = U22(0, 0);
	ux(2, 0) = dUdu11(0, 0); 	ux(2, 1) = dUdu12(0, 0);
	ux(3, 0) = dUdu21(0, 0); 	ux(3, 1) = dUdu22(0, 0);
	ux(0, 2) = dUdv11(0, 0);	ux(0, 3) = dUdv12(0, 0);
	ux(1, 2) = dUdv21(0, 0);	ux(1, 3) = dUdv22(0, 0);
	ux(2,2) = 0; ux(2,3) = 0; ux(3,2) = 0; ux(3,3) = 0;

	uy(0, 0) = U11(1, 0);		uy(0, 1) = U12(1, 0);
	uy(1, 0) = U21(1, 0);	 	uy(1, 1) = U22(1, 0);
	uy(2, 0) = dUdu11(1, 0); 	uy(2, 1) = dUdu12(1, 0);
	uy(3, 0) = dUdu21(1, 0); 	uy(3, 1) = dUdu22(1, 0);
	uy(0, 2) = dUdv11(1, 0);	uy(0, 3) = dUdv12(1, 0);
	uy(1, 2) = dUdv21(1, 0);	uy(1, 3) = dUdv22(1, 0);
	uy(2,2) = 0; uy(2,3) = 0; uy(3,2) = 0; uy(3,3) = 0;

	uz(0, 0) = U11(2, 0);		uz(0, 1) = U12(2, 0);
	uz(1, 0) = U21(2, 0);	 	uz(1, 1) = U22(2, 0);
	uz(2, 0) = dUdu11(2, 0); 	uz(2, 1) = dUdu12(2, 0);
	uz(3, 0) = dUdu21(2, 0); 	uz(3, 1) = dUdu22(2, 0);
	uz(0, 2) = dUdv11(2, 0);	uz(0, 3) = dUdv12(2, 0);
	uz(1, 2) = dUdv21(2, 0);	uz(1, 3) = dUdv22(2, 0);
	uz(2,2) = 0; uz(2,3) = 0; uz(3,2) = 0; uz(3,3) = 0;

	vnl_matrix_fixed<double, 4, 1> uu, uv;
	uu(0, 0) = h1(u);	uu(1, 0) = h2(u);	uu(2, 0) = h3(u);	uu(3, 0) = h4(u);
	uv(0, 0) = h1(v);	uv(1, 0) = h2(v);	uv(2, 0) = h3(v);	uv(3, 0) = h4(v);

	vnl_double_1x1 uxn = uu.transpose() * ux * uv;
	vnl_double_1x1 uyn = uu.transpose() * uy * uv;
	vnl_double_1x1 uzn = uu.transpose() * uz * uv;

	// (Su, Sv) = (r*(Uu Uv) - 2(pu pv)UtU)(I+UtU)^-1

	// Have S, compute Su, Sv, step, repeat

	vnl_double_3x1 S11, S12, S21, S22;
	vnl_double_3x1 B11, B12, B21, B22;

	vnl_double_4x4 rmat;

	rmat(0,0) = r11;	rmat(0,1) = r12;
	rmat(1,0) = r21;	rmat(1,1) = r22;
	rmat(2,0) = drdu11; rmat(2,1) = drdu12;
	rmat(3,0) = drdu21; rmat(3,1) = drdu22;
	rmat(0,2) = drdv11;	rmat(0,3) = drdv12;
	rmat(1,2) = drdv21; rmat(1,3) = drdv22;
	rmat(2,2) = 0; 		rmat(3,2) = 0;
	rmat(2,3) = 0; 		rmat(3,3) = 0;

	S11 = r11*U11;
	S21 = r21*U21;
	S12 = r12*U12;
	S22 = r22*U22;

	B11 = x11 + S11;
	B12 = x12 + S12;
	B21 = x21 + S21;
	B22 = x22 + S22;

	int numIter = 10;
	double curru = 0;
	double currv = 0;
	double du = u / numIter;
	double dv = v / numIter;

	vnl_matrix_fixed<double,4,1> humat, hvmat, hupmat, hvpmat;
	hupmat(0,0) = h1p(u);	hupmat(1,0) = h2p(u);	hupmat(2,0) = h3p(u);	hupmat(3,0) = h4p(u);
	hvpmat(0,0) = h1p(v);	hvpmat(1,0) = h2p(v);	hvpmat(2,0) = h3p(v);	hvpmat(3,0) = h4p(v);

	vnl_double_3x1 dUdu, dUdv;
	vnl_double_1x1 dUdux, dUduy, dUduz;
	vnl_double_1x1 dUdvx, dUdvy, dUdvz;

	vnl_double_3x1 dpdu, dpdv;

	vnl_matrix_fixed<double,2,3> Uu, Pu, Su, A;
	vnl_matrix_fixed<double,3,3> I, IUUt;
	I(0,0) = 1;	I(0,1) = 0;	I(0,2) = 0;
	I(1,0) = 0;	I(1,1) = 1; I(1,2) = 0;
	I(2,0) = 0; I(2,1) = 0; I(2,2) = 1;

	vnl_double_3x1 S, Sm, U;
	vnl_matrix_fixed<double,3,4> Smm;
	double r;

	double A00,A01,A02,A10,A11,A12,A20,A21,A22;
	double IU00,IU01,IU02,IU10,IU11,IU12,IU20,IU21,IU22;

	S.set_column(0,S11.get_column(0));
	U.set_column(0,U11.get_column(0));
	r = r11;

	Sm(0,0) = 0; Sm(1,0) = 0; Sm(2,0) = 0;

	vnl_double_3x1 Upu, Umu, Upv, Umv, Uuu, Uvv, temp1, temp2, temp3, temp4, tempU;

	for (int k = 0; k < 4; k++)
	{
		switch(k)
		{
		case 0:
			S.set_column(0,S11.get_column(0));
			U.set_column(0,U11.get_column(0));
			r = r11;

			curru = 0;
			currv = 0;
			du = u / numIter;
			dv = v / numIter;
			break;
		case 1:
			S.set_column(0,S12.get_column(0));
			U.set_column(0,U12.get_column(0));
			r = r12;

			curru = 0;
			currv = 1;
			du = u / numIter;
			dv = (v-1) / numIter;
			break;
		case 2:
			S.set_column(0,S21.get_column(0));
			U.set_column(0,U21.get_column(0));
			r = r21;

			curru = 1;
			currv = 0;
			du = (u-1) / numIter;
			dv = v / numIter;
			break;
		case 3:
			S.set_column(0,S22.get_column(0));
			U.set_column(0,U22.get_column(0));
			r = r22;

			curru = 1;
			currv = 1;
			du = (u-1) / numIter;
			dv = (v-1) / numIter;
			break;
		}


		double del = 1e-5; 

		for (int i = 0; i < numIter; i++)
		{
			humat(0,0) = h1(curru);	humat(1,0) = h2(curru);	humat(2,0) = h3(curru);	humat(3,0) = h4(curru);
			hvmat(0,0) = h1(currv);	hvmat(1,0) = h2(currv);	hvmat(2,0) = h3(currv);	hvmat(3,0) = h4(currv);

			hupmat(0,0) = h1p(curru); hupmat(1,0) = h2p(curru); hupmat(2,0) = h3p(curru); hupmat(3,0) = h4p(curru);
			hvpmat(0,0) = h1p(currv); hvpmat(1,0) = h2p(currv); hvpmat(2,0) = h3p(currv); hvpmat(3,0) = h4p(currv);

			if (onEdge)
			{

				bislerp(U11, U12, U21, U22, curru, currv, U);

				bislerp(U11, U12, U21, U22, curru+del, currv, Upu);
				bislerp(U11, U12, U21, U22, curru-del, currv, Umu);
				bislerp(U11, U12, U21, U22, curru, currv+del, Upv);
				bislerp(U11, U12, U21, U22, curru, currv-del, Umv);
			}
			else
			{
				// This should be replaced with actual closed-form computations
				temp1 = squad(U00,U01,U02,U03,curru+del);
				temp2 = squad(U10,U11,U12,U13,curru+del);
    			temp3 = squad(U20,U21,U22,U23,curru+del);
    			temp4 = squad(U30,U31,U32,U33,curru+del);
    			Upu = squad(temp1,temp2,temp3,temp4,currv);

    			temp1 = squad(U00,U01,U02,U03,curru-del);
				temp2 = squad(U10,U11,U12,U13,curru-del);
    			temp3 = squad(U20,U21,U22,U23,curru-del);
    			temp4 = squad(U30,U31,U32,U33,curru-del);
    			Umu = squad(temp1,temp2,temp3,temp4,currv);

    			temp1 = squad(U00,U01,U02,U03,curru);
				temp2 = squad(U10,U11,U12,U13,curru);
    			temp3 = squad(U20,U21,U22,U23,curru);
    			temp4 = squad(U30,U31,U32,U33,curru);
    			Upv = squad(temp1,temp2,temp3,temp4,currv+del);

    			temp1 = squad(U00,U01,U02,U03,curru);
				temp2 = squad(U10,U11,U12,U13,curru);
    			temp3 = squad(U20,U21,U22,U23,curru);
    			temp4 = squad(U30,U31,U32,U33,curru);
    			Umv = squad(temp1,temp2,temp3,temp4,currv-del);
			}


			Uuu = (Upu - Umu)/(2*del);
			Uvv = (Upv - Umv)/(2*del);

			Uu(0,0) = Uuu(0,0);
			Uu(0,1) = Uuu(1,0);
			Uu(0,2) = Uuu(2,0);
			Uu(1,0) = Uvv(0,0);
			Uu(1,1) = Uvv(1,0);
			Uu(1,2) = Uvv(2,0);

			Pu(0,0) = (hupmat.transpose() * hx * hvmat)(0,0);
			Pu(0,1) = (hupmat.transpose() * hy * hvmat)(0,0);
			Pu(0,2) = (hupmat.transpose() * hz * hvmat)(0,0);
			Pu(1,0) = (humat.transpose() * hx * hvpmat)(0,0);
			Pu(1,1) = (humat.transpose() * hy * hvpmat)(0,0);
			Pu(1,2) = (humat.transpose() * hz * hvpmat)(0,0);

			r = (humat.transpose() * rmat * hvmat)(0,0);

			IUUt = I + U*U.transpose();
//            IUUt = I - U*U.transpose();
			IU00 = IUUt(0,0); IU01 = IUUt(0,1); IU02 = IUUt(0,2);
			IU10 = IUUt(1,0); IU11 = IUUt(1,1); IU12 = IUUt(1,2);
			IU20 = IUUt(2,0); IU21 = IUUt(2,1); IU22 = IUUt(2,2);

			A = ((r*Uu)-(2.0*Pu*U*U.transpose()));
//            A = r*Uu;
//           Su = A-Pu*U*U.transpose();
			A00 = A(0,0); A01 = A(0,1); A02 = A(0,2);
			A10 = A(1,0); A11 = A(1,1); A12 = A(1,2);

			Su(0,0) = (A02*IU11*IU20 - A01*IU12*IU20 - A02*IU10*IU21 + A00*IU12*IU21 +
					A01*IU10*IU22 - A00*IU11*IU22)/(IU02*IU11*IU20 - IU01*IU12*IU20 -
							IU02*IU10*IU21 + IU00*IU12*IU21 + IU01*IU10*IU22 - IU00*IU11*IU22);

			Su(0,1) = -1*(A02*IU01*IU20 - A01*IU02*IU20 - A02*IU00*IU21 + A00*IU02*IU21 +
					A01*IU00*IU22 - A00*IU01*IU22)/(IU02*IU11*IU20 - IU01*IU12*IU20 -
							IU02*IU10*IU21 + IU00*IU12*IU21 + IU01*IU10*IU22 - IU00*IU11*IU22);

			Su(0,2) = -1*(A02*IU01*IU10 - A01*IU02*IU10 - A02*IU00*IU11 + A00*IU02*IU11 +
					A01*IU00*IU12 - A00*IU01*IU12)/(-IU02*IU11*IU20 + IU01*IU12*IU20 +
							IU02*IU10*IU21 - IU00*IU12*IU21 - IU01*IU10*IU22 + IU00*IU11*IU22);

			Su(1,0) = -1*(-1*A12*IU11*IU20 + A11*IU12*IU20 + A12*IU10*IU21 - A10*IU12*IU21 -
					A11*IU10*IU22 + A10*IU11*IU22)/(IU02*IU11*IU20 - IU01*IU12*IU20 -
							IU02*IU10*IU21 + IU00*IU12*IU21 + IU01*IU10*IU22 - IU00*IU11*IU22);

			Su(1,1) = -1*(A12*IU01*IU20 - A11*IU02*IU20 - A12*IU00*IU21 + A10*IU02*IU21 +
					A11*IU00*IU22 - A10*IU01*IU22)/(IU02*IU11*IU20 - IU01*IU12*IU20 -
							IU02*IU10*IU21 + IU00*IU12*IU21 + IU01*IU10*IU22 - IU00*IU11*IU22);

			Su(1,2) = -1*(A12*IU01*IU10 - A11*IU02*IU10 - A12*IU00*IU11 + A10*IU02*IU11 +
					A11*IU00*IU12 - A10*IU01*IU12)/(-IU02*IU11*IU20 + IU01*IU12*IU20 +
							IU02*IU10*IU21 - IU00*IU12*IU21 - IU01*IU10*IU22 + IU00*IU11*IU22);

			S.set_column(0,S.get_column(0) + du*Su.get_row(0) + dv*Su.get_row(1));

			curru += du;
			currv += dv;
		}
		Sm = S;
		Smm.set_column(k,S.get_column(0));
	}

	Sm.set_column(0,(1-u)*(1-v)*Smm.get_column(0) + (1-u)*(v)*Smm.get_column(1) +
			(u)*(1-v)*Smm.get_column(2) + (u)*(v)*Smm.get_column(3) );
	Sm.set_column(0,Smm.get_column(0));

	r = (hu.transpose() * rmat * hv)(0,0);

	U.set_column(0,Sm.get_column(0).normalize());

    Vector3D newU = Vector3D(U(0,0),U(1,0),U(2,0));

	Vector3D newX = Vector3D(xn(0,0),yn(0,0),zn(0,0));

	//M3DSpoke spoke = M3DSpoke(xn(0,0),yn(0,0), zn(0,0), testU(0,0), testU(1,0), testU(2,0), r);
	M3DSpoke spoke = M3DSpoke(xn(0,0),yn(0,0), zn(0,0), U(0,0), U(1,0), U(2,0), r);
	std::cout << "S: " << U(0,0) << " " << U(1,0) << " " << U(2,0) << std::endl;
	std::cout << "R: " << r << std::endl;

	return spoke;

}

M3DSpoke M3DInterpolater::interpolateSpoke( M3DFigure *figure, double u, double v, int side )
{
		int uBase = ( int ) floor ( u );
	    int vBase = ( int ) floor ( v );

	    M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );

	    if ( uBase == tempFigure->getRowCount() - 1 )
	    {
	        uBase = uBase - 1;
	    }
	    if ( vBase == tempFigure->getColumnCount() - 1 )
	    {
	        vBase = vBase - 1;
	    }
	    // Get four corner atoms

	    M3DPrimitive *atom11 = tempFigure->getPrimitivePtr ( uBase,vBase );
	    M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*> ( atom11 );

	    M3DPrimitive *atom21 = tempFigure->getPrimitivePtr ( uBase+1,vBase );
	    M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*> ( atom21 );

	    M3DPrimitive *atom12 = tempFigure->getPrimitivePtr ( uBase,vBase+1 );
	    M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*> ( atom12 );

	    M3DPrimitive *atom22 = tempFigure->getPrimitivePtr ( uBase+1,vBase+1 );
	    M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*> ( atom22 );

	    Vector3D x11 = quadAtom11->getX();
	    Vector3D x21 = quadAtom21->getX();
	    Vector3D x12 = quadAtom12->getX();
	    Vector3D x22 = quadAtom22->getX();

	    // Corner 11

	    Vector3D l1 = x21 - x11;
	    Vector3D l2 = x12 - x11;

	    double mag1 = l1.norm();
	    double mag2 = l2.norm();

	    double dot = (l1 * l2) / (mag1 * mag2);
	    double t11 = acos(dot);

	    // Corner 12

	    l1 = x22 - x12;
		l2 = x11 - x12;

		mag1 = l1.norm();
		mag2 = l2.norm();

		dot = (l1 * l2) / (mag1 * mag2);
		double t12 = acos(dot);

	    // Corner 21

		l1 = x22 - x21;
		l2 = x11 - x21;

		mag1 = l1.norm();
		mag2 = l2.norm();

		dot = (l1 * l2) / (mag1 * mag2);
		double t21 = acos(dot);

	    // Corner 22

		l1 = x21 - x22;
		l2 = x12 - x22;

		mag1 = l1.norm();
		mag2 = l2.norm();

		dot = (l1 * l2) / (mag1 * mag2);
		double t22 = acos(dot);

		M3DSpoke spoke;

        // angle > 160 degree
		if ((t11 > 2.8) || (t12 > 2.8) || (t21 > 2.8) || (t22 > 2.8))
		{
			std::cout << "u: " << u << ", v: " << v << ", using hermite interpolation" << std::endl;
			spoke = interpolateSpoke_hermite( figure, u, v, side );
		}
		else
		{
			std::cout << "u: " << u << ", v: " << v << ", using spoke interpolation" << std::endl;
			spoke = interpolateSpoke_new( figure, u, v, side );
            //spoke = interpolate(figure, u, v, side);
		}

		//interpolateSpoke_new( figure, u, v, side );

		return spoke;

		//return interpolateSpoke_hermite(figure, u, v, side);
}

//M3DSpoke* M3DInterpolater::interpolateSpoke_new( M3DFigure *figure, double u, double v, int side )
//{
//	int uBase = (int)floor(u);
//	int vBase = (int)floor(v);
//
//	u = u - uBase;
//	v = v - vBase;
//
//	M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );
//
//	if ( uBase == tempFigure->getRowCount() - 1 )
//	{
//		uBase = uBase - 1;
//		u = 1;
//	}
//
//	if ( vBase == tempFigure->getColumnCount() - 1 )
//	{
//		vBase = vBase - 1;
//		v = 1;
//	}
//
//
//}


M3DSpoke M3DInterpolater::interpolateSpoke_mean( M3DFigure *figure, double u, double v, int side )
{
	//std::cout << "interp mean" << std::endl;
    int uBase = ( int ) floor ( u );
    int vBase = ( int ) floor ( v );

    u = u - uBase;
    v = v - vBase;

    M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );

    if ( uBase == tempFigure->getRowCount() - 1 )
    {
        uBase = uBase - 1;
        u=1;
    }
    if ( vBase == tempFigure->getColumnCount() - 1 )
    {
        vBase = vBase - 1;
        v=1;
    }
    // Get four corner atoms

    M3DPrimitive *atom11 = tempFigure->getPrimitivePtr ( uBase,vBase );
    M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*> ( atom11 );

    M3DPrimitive *atom21 = tempFigure->getPrimitivePtr ( uBase+1,vBase );
    M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*> ( atom21 );

    M3DPrimitive *atom12 = tempFigure->getPrimitivePtr ( uBase,vBase+1 );
    M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*> ( atom12 );

    M3DPrimitive *atom22 = tempFigure->getPrimitivePtr ( uBase+1,vBase+1 );
    M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*> ( atom22 );

    Vector3D x11 = quadAtom11->getX();
    Vector3D x21 = quadAtom21->getX();
    Vector3D x12 = quadAtom12->getX();
    Vector3D x22 = quadAtom22->getX();

    //std::cout << "x11: " << x11.getX() << ", " << x11.getY() << ", " << x11.getZ() << std::endl;
    //std::cout << "x21: " << x21.getX() << ", " << x21.getY() << ", " << x21.getZ() << std::endl;
    //std::cout << "x12: " << x12.getX() << ", " << x12.getY() << ", " << x12.getZ() << std::endl;
    //std::cout << "x22: " << x22.getX() << ", " << x22.getY() << ", " << x22.getZ() << std::endl;

    double r11, r12, r21, r22, ru11, ru12, ru21, ru22, rv11, rv12, rv21, rv22;

    Vector3D U0_11, U0_12, U0_21, U0_22, B0_11, B0_12, B0_21, B0_22, S0_11, S0_12,S0_21,S0_22;

    if ( side == 0 )
    {
        U0_11 = quadAtom11->getU0();
        U0_21 = quadAtom21->getU0();
        U0_12 = quadAtom12->getU0();
        U0_22 = quadAtom22->getU0();

        B0_11 = quadAtom11->getX() + quadAtom11->getR0() * quadAtom11->getU0();
        B0_21 = quadAtom21->getX() + quadAtom21->getR0() * quadAtom21->getU0();
        B0_12 = quadAtom12->getX() + quadAtom12->getR0() * quadAtom12->getU0();
        B0_22 = quadAtom22->getX() + quadAtom22->getR0() * quadAtom22->getU0();

        S0_11 = U0_11 * quadAtom11->getR0();
        S0_21 = U0_21 * quadAtom21->getR0();
        S0_12 = U0_12 * quadAtom12->getR0();
        S0_22 = U0_22 * quadAtom22->getR0();

        r11 = quadAtom11->getR0();
        r21 = quadAtom21->getR0();
        r12 = quadAtom12->getR0();
        r22 = quadAtom22->getR0();

//		std::cout << "r11: " << quadAtom11->getR0() << std::endl;
//		std::cout << "r21: " << quadAtom21->getR0() << std::endl;
//		std::cout << "r12: " << quadAtom12->getR0() << std::endl;
//		std::cout << "r22: " << quadAtom22->getR0() << std::endl;
    }
    else
    {
        U0_11 = quadAtom11->getU1();
        U0_21 = quadAtom21->getU1();
        U0_12 = quadAtom12->getU1();
        U0_22 = quadAtom22->getU1();

        B0_11 = quadAtom11->getX() + quadAtom11->getR1() * quadAtom11->getU1();
        B0_21 = quadAtom21->getX() + quadAtom21->getR1() * quadAtom21->getU1();
        B0_12 = quadAtom12->getX() + quadAtom12->getR1() * quadAtom12->getU1();
        B0_22 = quadAtom22->getX() + quadAtom22->getR1() * quadAtom22->getU1();

        S0_11 = U0_11 * quadAtom11->getR1();
        S0_21 = U0_21 * quadAtom21->getR1();
        S0_12 = U0_12 * quadAtom12->getR1();
        S0_22 = U0_22 * quadAtom22->getR1();

        r11 = quadAtom11->getR1();
        r21 = quadAtom21->getR1();
        r12 = quadAtom12->getR1();
        r22 = quadAtom22->getR1();
    }

    //std::cout << "U0_11: " << U0_11.getX() << ", " << U0_11.getY() << ", " << U0_11.getZ() << std::endl;
    //std::cout << "U0_21: " << U0_21.getX() << ", " << U0_21.getY() << ", " << U0_21.getZ() << std::endl;
    //std::cout << "U0_12: " << U0_12.getX() << ", " << U0_12.getY() << ", " << U0_12.getZ() << std::endl;
    //std::cout << "U0_22: " << U0_22.getX() << ", " << U0_22.getY() << ", " << U0_22.getZ() << std::endl;

    Vector3D du11, du21, du12, du22, dv11, dv21, dv12, dv22;

    Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;
    Vector3D dB0du11, dB0dv11, dB0du21, dB0dv21, dB0du12, dB0dv12, dB0du22, dB0dv22;
    Vector3D dS0du11, dS0dv11, dS0du21, dS0dv21, dS0du12, dS0dv12, dS0du22, dS0dv22;

    dB0du11 = B0_21 - B0_11;
    dB0dv11 = B0_12 - B0_11;

    dB0du21 = B0_21 - B0_11;
    dB0dv21 = B0_22 - B0_21;

    dB0du12 = B0_22 - B0_12;
    dB0dv12 = B0_12 - B0_11;

    dB0du22 = B0_22 - B0_12;
    dB0dv22 = B0_22 - B0_21;

    // Get derivatives for each corner, starting at 11, then 21, 12, and 22

    // Corner 11
    // du

    if ( uBase == 0 )
    {
        du11 = x21 - x11;
        ru11 = r21 - r11;
        dU0du11 = U0_21 - U0_11;
        dS0du11 = S0_21 - S0_11;
        dB0du11 = B0_21 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase-1, vBase ) );
        du11 = ( x21 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            ru11 = ( r21 - tempAtom->getR0() ) / 2;
            dU0du11 = ( U0_21 - tempAtom->getU0() ) / 2;
            dS0du11 = ( S0_21 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dB0du11 = ( B0_21 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            ru11 = ( r21 - tempAtom->getR1() ) / 2;
            dU0du11 = ( U0_21 - tempAtom->getU1() ) / 2;
            dS0du11 = ( S0_21 - ( tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
            dB0du11 = ( B0_21 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }

        //dU0du11.print();
    }

    //dv

    if ( vBase == 0 )
    {
        dv11 = x12 - x11;
        rv11 = r12 - r11;
        dU0dv11 = U0_12 - U0_11;
        dS0dv11 = S0_12 - S0_11;
        dB0dv11 = B0_12 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase, vBase-1 ) );
        dv11 = ( x12 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            rv11 = ( r12 - tempAtom->getR0() ) / 2;
            dU0dv11 = ( U0_12 - tempAtom->getU0() ) / 2;
            dS0dv11 = ( U0_12 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dB0dv11 = ( B0_12 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            rv11 = ( r12 - tempAtom->getR1() ) / 2;
            dU0dv11 = ( U0_12 - tempAtom->getU1() ) / 2;
            dS0dv11 = ( U0_12 - ( tempAtom->getR0() * tempAtom->getU1() ) ) / 2;
            dB0dv11 = ( B0_12 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }
    }

    // Corner 21

    // du

    //std::cout << tempFigure->getRowCount() - 1 << std::endl;
    //std::cout << tempFigure->getColumnCount() - 1  << std::endl;
    if ( uBase + 1 == tempFigure->getRowCount() - 1 )
    {
        du21 = x21 - x11;
        ru21 = r21 - r11;
        dS0du21 = S0_21 - S0_11;
        dU0du21 = U0_21 - U0_11;
        dB0du21 = B0_21 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+2, vBase ) );
        du21 = ( tempAtom->getX() - x11 ) / 2;

        if ( side == 0 )
        {
            ru21 = ( tempAtom->getR0() - r11 ) / 2;
            dS0du21 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_11 ) / 2;
            dU0du21 = ( tempAtom->getU0() - U0_11 ) / 2;
            dB0du21 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_11 ) / 2;
        }
        else
        {
            ru21 = ( tempAtom->getR1() - r11 ) / 2;
            dS0du21 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_11 ) / 2;
            dU0du21 = ( tempAtom->getU1() - U0_11 ) / 2;
            dB0du21 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_11 ) / 2;
        }
    }

    //dv

    if ( vBase == 0 )
    {
        dv21 = x22 - x21;
        rv21 = r22 - r21;
        dS0dv21 = S0_22 - S0_21;
        dU0dv21 = U0_22 - U0_21;
        dB0dv21 = B0_22 - B0_21;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+1, vBase-1 ) );
        dv21 = ( x22 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            rv21 = ( r22 - tempAtom->getR0() ) / 2;
            dS0dv21 = ( S0_22 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dU0dv21 = ( U0_22 - tempAtom->getU0() ) / 2;
            dB0dv21 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            rv21 = ( r22 - tempAtom->getR1() ) / 2;
            dS0dv21 = ( S0_22 - ( tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
            dU0dv21 = ( U0_22 - tempAtom->getU1() ) / 2;
            dB0dv21 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }
    }

    // Corner 12

    // du

    if ( uBase == 0 )
    {
        du12 = x22 - x12;
        ru12 = r22 - r12;
        dS0du12 = S0_22 - S0_12;
        dU0du12 = U0_22 - U0_12;
        dB0du12 = B0_22 - B0_12;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase-1,vBase+1 ) );
        du12 = ( x22 - tempAtom->getX() ) / 2;

        if ( side == 0 )
        {
            ru12 = ( r22 - tempAtom->getR0() ) / 2;
            dS0du12 = ( S0_22 - ( tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
            dU0du12 = ( U0_22 - tempAtom->getU0() ) / 2;
            dB0du12 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) ) / 2;
        }
        else
        {
            ru12 = ( r22 - tempAtom->getR1() ) / 2;
            dS0du12 = ( S0_22 - ( tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
            dU0du12 = ( U0_22 - tempAtom->getU1() ) / 2;
            dB0du12 = ( B0_22 - ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) ) / 2;
        }
    }

    if ( vBase + 1 == tempFigure->getColumnCount() - 1 )
    {
        dv12 = x12 - x11;
        rv12 = r12 - r11;
        dS0dv12 = S0_12 - S0_11;
        dU0dv12 = U0_12 - U0_11;
        dB0dv12 = B0_12 - B0_11;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase, vBase+2 ) );
        dv12 = ( tempAtom->getX() - x11 ) / 2;

        if ( side == 0 )
        {
            rv12 = ( tempAtom->getR0() - r11 ) / 2;
            dS0dv12 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_11 ) / 2;
            dU0dv12 = ( tempAtom->getU0() - U0_11 ) / 2;
            dB0dv12 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_11 ) / 2;
        }
        else
        {
            rv12 = ( tempAtom->getR1() - r11 ) / 2;
            dS0dv12 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_11 ) / 2;
            dU0dv12 = ( tempAtom->getU1() - U0_11 ) / 2;
            dB0dv12 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_11 ) / 2;
        }
    }

    // Corner 22

    if ( uBase + 1 == tempFigure->getRowCount() - 1 )
    {
        du22 = x22 - x12;
        ru22 = r22 - r12;
        dS0du22 = S0_22 - S0_12;
        dU0du22 = U0_22 - U0_12;
        dB0du22 = B0_22 - B0_12;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+2,vBase+1 ) );
        du22 = ( tempAtom->getX() - x12 ) / 2;

        if ( side == 0 )
        {
            ru22 = ( tempAtom->getR0() - r12 ) / 2;
            dS0du22 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_12 ) / 2;
            dU0du22 = ( tempAtom->getU0() - U0_12 ) / 2;
            dB0du22 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_12 ) / 2;
        }
        else
        {
            ru22 = ( tempAtom->getR1() - r12 ) / 2;
            dS0du22 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_12 ) / 2;
            dU0du22 = ( tempAtom->getU1() - U0_12 ) / 2;
            dB0du22 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_12 ) / 2;
        }
    }

    if ( vBase + 1 == tempFigure->getColumnCount() - 1 )
    {
        dv22 = x22 - x21;
        rv22 = r22 - r21;
        dS0dv22 = S0_22 - S0_21;
        dU0dv22 = U0_22 - U0_21;
        dB0dv22 = B0_22 - B0_21;
    }
    else
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( uBase+1,vBase+2 ) );
        dv22 = ( tempAtom->getX() - x21 ) / 2;

        if ( side == 0 )
        {
            rv22 = ( tempAtom->getR0() - r22 ) / 2;
            dS0dv22 = ( ( tempAtom->getR0() * tempAtom->getU0() ) - S0_21 ) / 2;
            dU0dv22 = ( tempAtom->getU0() - U0_21 ) / 2;
            dB0dv22 = ( ( tempAtom->getX() + tempAtom->getR0() * tempAtom->getU0() ) - B0_21 ) / 2;
        }
        else
        {
            rv22 = ( tempAtom->getR1() - r22 ) / 2;
            dS0dv22 = ( ( tempAtom->getR1() * tempAtom->getU1() ) - S0_21 ) / 2;
            dU0dv22 = ( tempAtom->getU1() - U0_21 ) / 2;
            dB0dv22 = ( ( tempAtom->getX() + tempAtom->getR1() * tempAtom->getU1() ) - B0_21 ) / 2;
        }
    }

    if (side == 0)
    {
		Vector3D dU0du11_2 = tempFigure->computeDU0du(uBase, vBase);

		//std::cout << "Old: " << dU0du11 << std::endl;
		//std::cout << "New: " << dU0du11_2 << std::endl;
    }

//    du11 = x21 - x11;
//    dv11 = x12 - x11;
//    dU0du11 = U0_21 - U0_11;
//    dU0dv11 = U0_12 - U0_11;
//    dS0du11 = S0_21 - S0_11;
//    dS0dv11 = S0_12 - S0_11;
//    ru11 = r21 - r11;
//    rv11 = r12 - r11;
////
////
//    du21 = x21 - x11;
//    dv21 = x22 - x21;
//    dU0du21 = U0_21 - U0_11;
//    dU0dv21 = U0_22 - U0_21;
//    dS0du21 = S0_21 - S0_11;
//    dS0dv21 = S0_22 - S0_21;
//    ru21 = r21 - r11;
//    rv21 = r22 - r21;
////
//    du12 = x22 - x12;
//    dv12 = x12 - x11;
//    dU0du12 = U0_22 - U0_12;
//    dU0dv12 = U0_12 - U0_11;
//    dS0du12 = S0_22 - S0_12;
//    dS0dv12 = S0_12 - S0_11;
//    ru12 = r22 - r12;
//    rv12 = r12 - r11;
////
//    du22 = x22 - x12;
//    dv22 = x22 - x21;
//    dU0du22 = U0_22 - U0_12;
//    dU0dv22 = U0_22 - U0_21;
//    dS0du22 = S0_22 - S0_12;
//    dS0dv22 = S0_22 - S0_21;
//    ru22 = r22 - r12;
//    rv22 = r22 - r21;



//	Vector3D du11 = (x21 - x01) / 2;
//	Vector3D du21 = (x31 - x11) / 2;
//	Vector3D du12 = (x22 - x02) / 2;
//	Vector3D du22 = (x32 - x12) / 2;
//
//	Vector3D dv11 = (x12 - x10) / 2;
//	Vector3D dv21 = (x22 - x20) / 2;
//	Vector3D dv12 = (x13 - x11) / 2;
//	Vector3D dv22 = (x23 - x21) / 2;

//	std::cout << "dU0du11: " << dU0du11.getX() << ", " << dU0du11.getY() << ", " << dU0du11.getZ() << std::endl;
//	std::cout << "dU0du21: " << dU0du21.getX() << ", " << dU0du21.getY() << ", " << dU0du21.getZ() << std::endl;
//	std::cout << "dU0du12: " << dU0du12.getX() << ", " << dU0du12.getY() << ", " << dU0du12.getZ() << std::endl;
//	std::cout << "dU0du22: " << dU0du22.getX() << ", " << dU0du22.getY() << ", " << dU0du22.getZ() << std::endl;
//
//	std::cout << "dU0dv11: " << dU0dv11.getX() << ", " << dU0dv11.getY() << ", " << dU0dv11.getZ() << std::endl;
//	std::cout << "dU0dv21: " << dU0dv21.getX() << ", " << dU0dv21.getY() << ", " << dU0dv21.getZ() << std::endl;
//	std::cout << "dU0dv12: " << dU0dv12.getX() << ", " << dU0dv12.getY() << ", " << dU0dv12.getZ() << std::endl;
//	std::cout << "dU0dv22: " << dU0dv22.getX() << ", " << dU0dv22.getY() << ", " << dU0dv22.getZ() << std::endl;

    // Get unit normals for the four corners and project derivatives on to the tangent plane
    Vector3D n11 = quadAtom11->getU0() - quadAtom11->getU1();
    n11.normalize();
    Vector3D n21 = quadAtom21->getU0() - quadAtom21->getU1();
    n21.normalize();
    Vector3D n12 = quadAtom12->getU0() - quadAtom12->getU1();
    n12.normalize();
    Vector3D n22 = quadAtom22->getU0() - quadAtom22->getU1();
    n22.normalize();

//	n11 = du11.cross(dv11);
//	n11.normalize();
//	n12 = du12.cross(dv12);
//	n12.normalize();
//	n21 = du21.cross(dv21);
//	n21.normalize();
//	n22 = du22.cross(dv22);
//	n22.normalize();

    Vector3D du11t = du11 - ( du11 * n11 ) * n11;
    Vector3D du21t = du21 - ( du21 * n21 ) * n21;
    Vector3D du12t = du12 - ( du12 * n12 ) * n12;
    Vector3D du22t = du22 - ( du22 * n22 ) * n22;

    Vector3D dv11t = dv11 - ( dv11 * n11 ) * n11;
    Vector3D dv21t = dv21 - ( dv21 * n21 ) * n21;
    Vector3D dv12t = dv12 - ( dv12 * n12 ) * n12;
    Vector3D dv22t = dv22 - ( dv22 * n22 ) * n22;

    // Build matrices for hermite interpolation of medial sheet
    double hx[16] = { x11.getX(), x21.getX(), du11t.getX(), du21t.getX(),
                      x12.getX(), x22.getX(), du12t.getX(), du22t.getX(), dv11t.getX(),
                      dv21t.getX(), 0, 0, dv12t.getX(), dv22t.getX(), 0, 0
                    };

    double hy[16] = { x11.getY(), x21.getY(), du11t.getY(), du21t.getY(),
                      x12.getY(), x22.getY(), du12t.getY(), du22t.getY(), dv11t.getY(),
                      dv21t.getY(), 0, 0, dv12t.getY(), dv22t.getY(), 0, 0
                    };

    double hz[16] = { x11.getZ(), x21.getZ(), du11t.getZ(), du21t.getZ(),
                      x12.getZ(), x22.getZ(),du12t.getZ(), du22t.getZ(), dv11t.getZ(),
                      dv21t.getZ(), 0, 0, dv12t.getZ(), dv22t.getZ(), 0, 0
                    };

    double r_hermite[16] = { r11, r21, ru11, ru21, r12, r22, ru12, ru22, rv11, rv21, 0, 0, rv12, rv22, 0, 0};

    Matrix r_hermite_mat = Matrix ( 4,4,r_hermite, true );

    Matrix hxmat = Matrix ( 4, 4, hx, true );
    Matrix hymat = Matrix ( 4, 4, hy, true );
    Matrix hzmat = Matrix ( 4, 4, hz, true );

    //U0_11.print();

    double hu[4] = { h1 ( u ), h2 ( u ), h3 ( u ), h4 ( u ) };
    double hv[4] = { h1 ( v ), h2 ( v ), h3 ( v ), h4 ( v ) };
    Matrix humat = Matrix ( 1, 4, hu, true );
    Matrix hvmat = Matrix ( 4, 1, hv, true );

    Matrix xn = humat * hxmat * hvmat;
    Matrix yn = humat * hymat * hvmat;
    Matrix zn = humat * hzmat * hvmat;

    double Bx[16] = { B0_11.getX(), B0_21.getX(), dB0du11.getX(), dB0du21.getX(),
                      B0_12.getX(), B0_22.getX(), dB0du12.getX(), dB0du22.getX(), dB0dv11.getX(),
                      dB0dv21.getX(), 0, 0, dB0dv12.getX(), dB0dv22.getX(), 0, 0
                    };

    double By[16] = { B0_11.getY(), B0_21.getY(), dB0du11.getY(), dB0du21.getY(),
                      B0_12.getY(), B0_22.getY(), dB0du12.getY(), dB0du22.getY(), dB0dv11.getY(),
                      dB0dv21.getY(), 0, 0, dB0dv12.getY(), dB0dv22.getY(), 0, 0
                    };

    double Bz[16] = { B0_11.getZ(), B0_21.getZ(), dB0du11.getZ(), dB0du21.getZ(),
                      B0_12.getZ(), B0_22.getZ(), dB0du12.getZ(), dB0du22.getZ(), dB0dv11.getZ(),
                      dB0dv21.getZ(), 0, 0, dB0dv12.getZ(), dB0dv22.getZ(), 0, 0
                    };

    Matrix bxmat = Matrix ( 4, 4, Bx, true );
    Matrix bymat = Matrix ( 4, 4, By, true );
    Matrix bzmat = Matrix ( 4, 4, Bz, true );

    Matrix Bxn = humat * bxmat * hvmat;
    Matrix Byn = humat * bymat * hvmat;
    Matrix Bzn = humat * bzmat * hvmat;


    // Calculate sRad matrices using finite differences
    /*Vector3D dU0du11 = (quadAtom21->getU0() - quadAtom01->getU0()) / 2;
    Vector3D dU0dv11 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

    Vector3D dU0du21 = (quadAtom31->getU0() - quadAtom11->getU0()) / 2;
    Vector3D dU0dv21 = (quadAtom22->getU0() - quadAtom20->getU0()) / 2;

    Vector3D dU0du12 = (quadAtom22->getU0() - quadAtom02->getU0()) / 2;
    Vector3D dU0dv12 = (quadAtom12->getU0() - quadAtom10->getU0()) / 2;

    Vector3D dU0du22 = (quadAtom32->getU0() - quadAtom12->getU0()) / 2;
    Vector3D dU0dv22 = (quadAtom23->getU0() - quadAtom21->getU0()) / 2;*/

    /*Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;

    if (side == 0)
    {
    	dU0du11 = quadAtom21->getU0() - quadAtom11->getU0();
    	dU0dv11 = quadAtom12->getU0() - quadAtom11->getU0();

    	dU0du21 = quadAtom21->getU0() - quadAtom11->getU0();
    	dU0dv21 = quadAtom22->getU0() - quadAtom21->getU0();

    	dU0du12 = quadAtom22->getU0() - quadAtom12->getU0();
    	dU0dv12 = quadAtom12->getU0() - quadAtom11->getU0();

    	dU0du22 = quadAtom22->getU0() - quadAtom12->getU0();
    	dU0dv22 = quadAtom22->getU0() - quadAtom21->getU0();
    }
    else
    {
    	dU0du11 = quadAtom21->getU1() - quadAtom11->getU1();
    	dU0dv11 = quadAtom12->getU1() - quadAtom11->getU1();

    	dU0du21 = quadAtom21->getU1() - quadAtom11->getU1();
    	dU0dv21 = quadAtom22->getU1() - quadAtom21->getU1();

    	dU0du12 = quadAtom22->getU1() - quadAtom12->getU1();
    	dU0dv12 = quadAtom12->getU1() - quadAtom11->getU1();

    	dU0du22 = quadAtom22->getU1() - quadAtom12->getU1();
    	dU0dv22 = quadAtom22->getU1() - quadAtom21->getU1();
    }*/

    // These form the non-orthogonal medial coordinate system for the spoke derivatives

    //Vector3D U0_11, U0_12, U0_21, U0_22;

    //if (side == 0)
    //{
    //	U0_11 = quadAtom11->getU0();
    //	U0_21 = quadAtom21->getU0();
    //	U0_12 = quadAtom12->getU0();
    //	U0_22 = quadAtom22->getU0();
    //}
    //else
    //{
    //	U0_11 = quadAtom11->getU1();
    //	U0_21 = quadAtom21->getU1();
    //	U0_12 = quadAtom12->getU1();
    //	U0_22 = quadAtom22->getU1();
    //}
    Vector3D du11p = n11.cross ( du11 );
    Vector3D du12p = n12.cross ( du12 );
    Vector3D du21p = n21.cross ( du21 );
    Vector3D du22p = n22.cross ( du22 );

    Vector3D du11tp = n11.cross ( du11t );
    Vector3D du12tp = n12.cross ( du12t );
    Vector3D du21tp = n21.cross ( du21t );
    Vector3D du22tp = n22.cross ( du22t );

    Vector3D du11t_norm = du11t / du11t.norm();
    Vector3D du12t_norm = du12t / du12t.norm();
    Vector3D du21t_norm = du21t / du21t.norm();
    Vector3D du22t_norm = du22t / du22t.norm();

    Vector3D du11pt_norm = du11tp / du11tp.norm();
    Vector3D du12pt_norm = du12tp / du12tp.norm();
    Vector3D du21pt_norm = du21tp / du21tp.norm();
    Vector3D du22pt_norm = du22tp / du22tp.norm();

    double basis11[9] = { du11t.getX(), du11t.getY(), du11t.getZ(), dv11t.getX(), dv11t.getY(), dv11t.getZ(),
                          U0_11.getX(), U0_11.getY(), U0_11.getZ()
                        };

    double basis21[9] = { du21t.getX(), du21t.getY(), du21t.getZ(), dv21t.getX(), dv21t.getY(), dv21t.getZ(),
                          U0_21.getX(), U0_21.getY(), U0_21.getZ()
                        };

    double basis12[9] = { du12t.getX(), du12t.getY(), du12t.getZ(), dv12t.getX(), dv12t.getY(), dv12t.getZ(),
                          U0_12.getX(), U0_12.getY(), U0_12.getZ()
                        };

    double basis22[9] = { du22t.getX(), du22t.getY(), du22t.getZ(), dv22t.getX(), dv22t.getY(), dv22t.getZ(),
                          U0_22.getX(), U0_22.getY(), U0_22.getZ()
                        };

    /*double basis11[9] = { du11t_norm.getX(), du11t_norm.getY(), du11t_norm.getZ(), du11pt_norm.getX(), du11pt_norm.getY(), du11pt_norm.getZ(),
    	U0_11.getX(), U0_11.getY(), U0_11.getZ() };

    double basis21[9] = { du21t_norm.getX(), du21t_norm.getY(), du21t_norm.getZ(), du21pt_norm.getX(), du21pt_norm.getY(), du21pt_norm.getZ(),
    	U0_21.getX(), U0_21.getY(), U0_21.getZ() };

    double basis12[9] = { du12t_norm.getX(), du12t_norm.getY(), du12t_norm.getZ(), du12pt_norm.getX(), du12pt_norm.getY(), du12pt_norm.getZ(),
    	U0_12.getX(), U0_12.getY(), U0_12.getZ() };

    double basis22[9] = { du22t_norm.getX(), du22t_norm.getY(), du22t_norm.getZ(), du22pt_norm.getX(), du22pt_norm.getY(), du22pt_norm.getZ(),
    	U0_22.getX(), U0_22.getY(), U0_22.getZ() };*/



    Matrix C11 = Matrix ( 3, 3, basis11, true );
    Matrix C21 = Matrix ( 3, 3, basis21, true );
    Matrix C12 = Matrix ( 3, 3, basis12, true );
    Matrix C22 = Matrix ( 3, 3, basis22, true );

    Matrix Ci11, Ci21, Ci12, Ci22;

    C11.inverse ( Ci11 );
    C21.inverse ( Ci21 );
    C12.inverse ( Ci12 );
    C22.inverse ( Ci22 );

    double dU0du11_coeffs[3] = { dU0du11.getX(), dU0du11.getY(), dU0du11.getZ() };
    double dU0dv11_coeffs[3] = { dU0dv11.getX(), dU0dv11.getY(), dU0dv11.getZ() };
    double dU0du21_coeffs[3] = { dU0du21.getX(), dU0du21.getY(), dU0du21.getZ() };
    double dU0dv21_coeffs[3] = { dU0dv21.getX(), dU0dv21.getY(), dU0dv21.getZ() };
    double dU0du12_coeffs[3] = { dU0du12.getX(), dU0du12.getY(), dU0du12.getZ() };
    double dU0dv12_coeffs[3] = { dU0dv12.getX(), dU0dv12.getY(), dU0dv12.getZ() };
    double dU0du22_coeffs[3] = { dU0du22.getX(), dU0du22.getY(), dU0du22.getZ() };
    double dU0dv22_coeffs[3] = { dU0dv22.getX(), dU0dv22.getY(), dU0dv22.getZ() };

    double dS0du11_coeffs[3] = { dS0du11.getX(), dS0du11.getY(), dS0du11.getZ() };
    double dS0dv11_coeffs[3] = { dS0dv11.getX(), dS0dv11.getY(), dS0dv11.getZ() };
    double dS0du21_coeffs[3] = { dS0du21.getX(), dS0du21.getY(), dS0du21.getZ() };
    double dS0dv21_coeffs[3] = { dS0dv21.getX(), dS0dv21.getY(), dS0dv21.getZ() };
    double dS0du12_coeffs[3] = { dS0du12.getX(), dS0du12.getY(), dS0du12.getZ() };
    double dS0dv12_coeffs[3] = { dS0dv12.getX(), dS0dv12.getY(), dS0dv12.getZ() };
    double dS0du22_coeffs[3] = { dS0du22.getX(), dS0du22.getY(), dS0du22.getZ() };
    double dS0dv22_coeffs[3] = { dS0dv22.getX(), dS0dv22.getY(), dS0dv22.getZ() };

    Matrix Au11 = Matrix ( 3,1,dU0du11_coeffs,true );
    Matrix Av11 = Matrix ( 3,1,dU0dv11_coeffs,true );
    Matrix Au21 = Matrix ( 3,1,dU0du21_coeffs,true );
    Matrix Av21 = Matrix ( 3,1,dU0dv21_coeffs,true );
    Matrix Au12 = Matrix ( 3,1,dU0du12_coeffs,true );
    Matrix Av12 = Matrix ( 3,1,dU0dv12_coeffs,true );
    Matrix Au22 = Matrix ( 3,1,dU0du22_coeffs,true );
    Matrix Av22 = Matrix ( 3,1,dU0dv22_coeffs,true );

    Matrix SAu11 = Matrix ( 3,1,dS0du11_coeffs,true );
    Matrix SAv11 = Matrix ( 3,1,dS0dv11_coeffs,true );
    Matrix SAu21 = Matrix ( 3,1,dS0du21_coeffs,true );
    Matrix SAv21 = Matrix ( 3,1,dS0dv21_coeffs,true );
    Matrix SAu12 = Matrix ( 3,1,dS0du12_coeffs,true );
    Matrix SAv12 = Matrix ( 3,1,dS0dv12_coeffs,true );
    Matrix SAu22 = Matrix ( 3,1,dS0du22_coeffs,true );
    Matrix SAv22 = Matrix ( 3,1,dS0dv22_coeffs,true );

    Matrix Bu11 = Ci11 * Au11;
    Matrix Bv11 = Ci11 * Av11;
    Matrix Bu21 = Ci21 * Au21;
    Matrix Bv21 = Ci21 * Av21;
    Matrix Bu12 = Ci12 * Au12;
    Matrix Bv12 = Ci12 * Av12;
    Matrix Bu22 = Ci22 * Au22;
    Matrix Bv22 = Ci22 * Av22;

    Matrix SBu11 = Ci11 * SAu11;
    Matrix SBv11 = Ci11 * SAv11;
    Matrix SBu21 = Ci21 * SAu21;
    Matrix SBv21 = Ci21 * SAv21;
    Matrix SBu12 = Ci12 * SAu12;
    Matrix SBv12 = Ci12 * SAv12;
    Matrix SBu22 = Ci22 * SAu22;
    Matrix SBv22 = Ci22 * SAv22;

    double b11[4] = { -1*Bu11 ( 0,0 ), -1*Bu11 ( 1,0 ), -1*Bv11 ( 0,0 ), -1*Bv11 ( 1,0 ) };
    double b21[4] = { -1*Bu21 ( 0,0 ), -1*Bu21 ( 1,0 ), -1*Bv21 ( 0,0 ), -1*Bv21 ( 1,0 ) };
    double b12[4] = { -1*Bu12 ( 0,0 ), -1*Bu12 ( 1,0 ), -1*Bv12 ( 0,0 ), -1*Bv12 ( 1,0 ) };
    double b22[4] = { -1*Bu22 ( 0,0 ), -1*Bu22 ( 1,0 ), -1*Bv22 ( 0,0 ), -1*Bv22 ( 1,0 ) };

    double Sb11[4] = { -1*SBu11 ( 0,0 ), -1*SBu11 ( 1,0 ), -1*SBv11 ( 0,0 ), -1*SBv11 ( 1,0 ) };
    double Sb21[4] = { -1*SBu21 ( 0,0 ), -1*SBu21 ( 1,0 ), -1*SBv21 ( 0,0 ), -1*SBv21 ( 1,0 ) };
    double Sb12[4] = { -1*SBu12 ( 0,0 ), -1*SBu12 ( 1,0 ), -1*SBv12 ( 0,0 ), -1*SBv12 ( 1,0 ) };
    double Sb22[4] = { -1*SBu22 ( 0,0 ), -1*SBu22 ( 1,0 ), -1*SBv22 ( 0,0 ), -1*SBv22 ( 1,0 ) };

    Matrix Srad11 = Matrix ( 2,2,b11,true );
    Matrix Srad21 = Matrix ( 2,2,b21,true );
    Matrix Srad12 = Matrix ( 2,2,b12,true );
    Matrix Srad22 = Matrix ( 2,2,b22,true );

    Matrix SSrad11 = Matrix ( 2,2,Sb11,true );
    Matrix SSrad21 = Matrix ( 2,2,Sb21,true );
    Matrix SSrad12 = Matrix ( 2,2,Sb12,true );
    Matrix SSrad22 = Matrix ( 2,2,Sb22,true );

//	std::cout << "SSrad: " << std::endl;
//	SSrad12.print();

    Matrix rSrad11 = Srad11 * quadAtom11->getR0();
    Matrix rSrad21 = Srad21 * quadAtom21->getR0();
    Matrix rSrad12 = Srad12 * quadAtom12->getR0();
    Matrix rSrad22 = Srad22 * quadAtom22->getR0();

//	std::cout << "rSrad: " << std::endl;
//	rSrad12.print();

    double Idat[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    Matrix Imat = Matrix ( 3,3,Idat,true );

    //11

    double pdat11[6] = {du11t.getX(), dv11t.getX(), du11t.getY(), dv11t.getY(), du11t.getZ(), dv11t.getZ() };
//	double pdat11[6] = {du11t.getX(), du11tp.getX(), du11t.getY(), du11tp.getY(), du11t.getZ(), du11tp.getZ() };
    Matrix pmat11 = Matrix ( 2,3,pdat11,true );

    double Udat11[3] = {U0_11.getX(), U0_11.getY(), U0_11.getZ() };
    Matrix Umat11 = Matrix ( 1,3,Udat11,true );

    Matrix Q11 = pmat11 * ( ( Umat11.t() * Umat11 ) - Imat );

    Matrix QQ11 = Q11 * Q11.t();
    Matrix QQi11;
    QQ11.inverse ( QQi11 );

    double rdat11[2] = {ru11, rv11};
    Matrix rmat11 = Matrix ( 2,1,rdat11,true );

    Matrix r2mat11 = pmat11 * Umat11.t();

    double dSdu11dat[6] = { dS0du11.getX(), dS0dv11.getX(), dS0du11.getY(), dS0dv11.getY(), dS0du11.getZ(), dS0dv11.getZ() };
    Matrix dSdu11mat = Matrix ( 2,3,dSdu11dat,true );

    double dUdu11dat[6] = { dU0du11.getX(), dU0dv11.getX(), dU0du11.getY(), dU0dv11.getY(), dU0du11.getZ(), dU0dv11.getZ() };
    Matrix dUdu11mat = Matrix ( 2,3,dUdu11dat,true );

    Matrix rS11 = ( dSdu11mat - ( rmat11*Umat11 ) ) * ( Q11.t() *QQi11 );

    Matrix testmult = Umat11 * Umat11.t();
//	testmult.print();

    //21

    double pdat21[6] = {du21t.getX(), dv21t.getX(), du21t.getY(), dv21t.getY(), du21t.getZ(), dv21t.getZ() };
//	double pdat21[6] = {du21t.getX(), du21tp.getX(), du21t.getY(), du21tp.getY(), du21t.getZ(), du21tp.getZ() };
    Matrix pmat21 = Matrix ( 2,3,pdat21,true );

    double Udat21[3] = {U0_21.getX(), U0_21.getY(), U0_21.getZ() };
    Matrix Umat21 = Matrix ( 1,3,Udat21,true );

    Matrix Q21 = pmat21 * ( ( Umat21.t() * Umat21 ) - Imat );

    Matrix QQ21 = Q21 * Q21.t();
    Matrix QQi21;
    QQ21.inverse ( QQi21 );

    double rdat21[2] = {ru21, rv21};
    Matrix rmat21 = Matrix ( 2,1,rdat21,true );

    Matrix r2mat21 = -1 * pmat21 * Umat21.t();

    double dSdu21dat[6] = { dS0du21.getX(), dS0dv21.getX(), dS0du21.getY(), dS0dv21.getY(), dS0du21.getZ(), dS0dv21.getZ() };
    Matrix dSdu21mat = Matrix ( 2,3,dSdu21dat,true );

    Matrix rS21 = ( dSdu21mat - ( rmat21*Umat21 ) ) * ( Q21.t() *QQi21 );

    //12

    double pdat12[6] = {du12t.getX(), dv12t.getX(), du12t.getY(), dv12t.getY(), du12t.getZ(), dv12t.getZ() };
//	double pdat12[6] = {du12t.getX(), du12tp.getX(), du12t.getY(), du12tp.getY(), du12t.getZ(), du12tp.getZ() };
    Matrix pmat12 = Matrix ( 2,3,pdat12,true );
    double Udat12[3] = {U0_12.getX(), U0_12.getY(), U0_12.getZ() };
    Matrix Umat12 = Matrix ( 1,3,Udat12,true );

    Matrix Q12 = pmat12 * ( ( Umat12.t() * Umat12 ) - Imat );

    Matrix QQ12 = Q12 * Q12.t();
    Matrix QQi12;
    QQ12.inverse ( QQi12 );

    double rdat12[2] = {ru12, rv12};
    Matrix rmat12 = Matrix ( 2,1,rdat12,true );

    Matrix r2mat12 = -1 * pmat12 * Umat12.t();

    double dSdu12dat[6] = { dS0du12.getX(), dS0dv12.getX(), dS0du12.getY(), dS0dv12.getY(), dS0du12.getZ(), dS0dv12.getZ() };
    Matrix dSdu12mat = Matrix ( 2,3,dSdu12dat,true );

    Matrix rS12 = ( dSdu12mat - ( rmat12*Umat12 ) ) * ( Q12.t() *QQi12 );

    //22

    double pdat22[6] = {du22t.getX(), dv22t.getX(), du22t.getY(), dv22t.getY(), du22t.getZ(), dv22t.getZ() };
//	double pdat22[6] = {du22t.getX(), du22tp.getX(), du22t.getY(), du22tp.getY(), du22t.getZ(), du22tp.getZ() };
    Matrix pmat22 = Matrix ( 2,3,pdat22,true );

    double Udat22[3] = {U0_22.getX(), U0_22.getY(), U0_22.getZ() };
    Matrix Umat22 = Matrix ( 1,3,Udat22,true );

    Matrix Q22 = pmat22 * ( ( Umat22.t() * Umat22 ) - Imat );

    Matrix QQ22 = Q22 * Q22.t();
    Matrix QQi22;
    QQ22.inverse ( QQi22 );

    double rdat22[2] = {ru22, rv22};
    Matrix rmat22 = Matrix ( 2,1,rdat22,true );

    Matrix r2mat22 = -1 * pmat22 * Umat22.t();

    double dSdu22dat[6] = { dS0du22.getX(), dS0dv22.getX(), dS0du22.getY(), dS0dv22.getY(), dS0du22.getZ(), dS0dv22.getZ() };
    Matrix dSdu22mat = Matrix ( 2,3,dSdu22dat,true );

    Matrix rS22 = ( dSdu22mat - ( rmat22*Umat22 ) ) * ( Q22.t() *QQi22 );

    //rSrad11 = SSrad11;
    //rSrad21 = SSrad21;
    //rSrad12 = SSrad12;
    //rSrad22 = SSrad22;

    //std::cout << "srads" << std::endl;
    //rSrad11.print();
    //rSrad21.print();
    //rSrad12.print();
    //rSrad22.print();

    Matrix rS12_2 = ( dSdu12mat + ( r2mat12*Umat12 ) ) * ( Q12.t() *QQi12 );

//	std::cout << "rSrad from formula, with pu/pv: " << std::endl;
//	rS12_2.t().print();
//
//	std::cout << "rSrad from formula, with ru/rv: " << std::endl;
//	rS12.t().print();

//	rmat11.print();
//	r2mat11.print();

    rSrad11 = rS11.t();
    rSrad12 = rS12.t();
    rSrad21 = rS21.t();
    rSrad22 = rS22.t();


    Matrix rSrad11RT = rSrad11 * rSrad11.t();
    Matrix rSrad21RT = rSrad21 * rSrad21.t();
    Matrix rSrad12RT = rSrad12 * rSrad12.t();
    Matrix rSrad22RT = rSrad22 * rSrad22.t();

    Matrix rSrad11TR = rSrad11.t() * rSrad11;
    Matrix rSrad21TR = rSrad21.t() * rSrad21;
    Matrix rSrad12TR = rSrad12.t() * rSrad12;
    Matrix rSrad22TR = rSrad22.t() * rSrad22;

//	cout << endl << "rSrad: " << endl;
//	cout << rSrad11(0,0) << ", " << rSrad11(0,1) << ", " << rSrad11(1,0) << ", " << rSrad11(1,1) << endl;
//	cout << rSrad21(0,0) << ", " << rSrad21(0,1) << ", " << rSrad21(1,0) << ", " << rSrad21(1,1) << endl;
//	cout << rSrad12(0,0) << ", " << rSrad12(0,1) << ", " << rSrad12(1,0) << ", " << rSrad12(1,1) << endl;
//	cout << rSrad22(0,0) << ", " << rSrad22(0,1) << ", " << rSrad22(1,0) << ", " << rSrad22(1,1) << endl;

    Vector L11, L21, L12, L22;
    Matrix V11, V21, V12, V22;

    Vector L1_11, L2_11, L1_21, L2_21, L1_12, L2_12, L1_22, L2_22;
    Matrix V1_11, V2_11, V1_21, V2_21, V1_12, V2_12, V1_22, V2_22;

    rSrad11.factorEV ( L11, V11, NON_SYM );
    rSrad21.factorEV ( L21, V21, NON_SYM );
    rSrad12.factorEV ( L12, V12, NON_SYM );
    rSrad22.factorEV ( L22, V22, NON_SYM );

    rSrad11RT.factorEV ( L1_11, V1_11, GENERAL );
    rSrad21RT.factorEV ( L1_21, V1_21, GENERAL );
    rSrad12RT.factorEV ( L1_12, V1_12, GENERAL );
    rSrad22RT.factorEV ( L1_22, V1_22, GENERAL );

    rSrad11TR.factorEV ( L2_11, V2_11, GENERAL );
    rSrad21TR.factorEV ( L2_21, V2_21, GENERAL );
    rSrad12TR.factorEV ( L2_12, V2_12, GENERAL );
    rSrad22TR.factorEV ( L2_22, V2_22, GENERAL );

//	std::cout << "eigvals" << std::endl;
//
//	L1_11.print();
//	L2_11.print();
//
//	L1_12.print();
//	L2_12.print();
//
//	L1_21.print();
//	L2_21.print();
//
//	L1_22.print();
//	L2_22.print();

    // Process corners to find out true eigenvalues

    double d1, d2, nd1, nd2, res1, res2, res3, res4, minres;
    Matrix dPP, dPN, dNP, dNN, test1, test2, test3, test4, temp1, temp2, temp3, temp4, EVs_11, EVs_21, EVs_12, EVs_22;

    //Corner 11

    d1 = sqrt ( L2_11 ( 0 ) );
    d2 = sqrt ( L2_11 ( 1 ) );

    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_11 * dPP * V2_11;
    test2 = V1_11 * dPN * V2_11;
    test3 = V1_11 * dNP * V2_11;
    test4 = V1_11 * dNN * V2_11;

    temp1 = test1 - rSrad11;
    temp2 = test2 - rSrad11;
    temp3 = test3 - rSrad11;
    temp4 = test4 - rSrad11;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_11 = dPP;

    if ( res2 < minres )
    {
        minres = res2;
        EVs_11 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_11 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_11 = dNN;
    }

    //std::cout << minres << std::endl;
    //Corner 21

    d1 = sqrt ( L2_21 ( 0 ) );
    d2 = sqrt ( L2_21 ( 1 ) );
    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_21 * dPP * V2_21;
    test2 = V1_21 * dPN * V2_21;
    test3 = V1_21 * dNP * V2_21;
    test4 = V1_21 * dNN * V2_21;

    temp1 = test1 - rSrad21;
    temp2 = test2 - rSrad21;
    temp3 = test3 - rSrad21;
    temp4 = test4 - rSrad21;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_21 = dPP;

    if ( res2 < minres )
    {
        minres = res2;
        EVs_21 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_21 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_21 = dNN;
    }

    //Corner 12

    d1 = sqrt ( L2_12 ( 0 ) );
    d2 = sqrt ( L2_12 ( 1 ) );
    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_12 * dPP * V2_12;
    test2 = V1_12 * dPN * V2_12;
    test3 = V1_12 * dNP * V2_12;
    test4 = V1_12 * dNN * V2_12;

    temp1 = test1 - rSrad12;
    temp2 = test2 - rSrad12;
    temp3 = test3 - rSrad12;
    temp4 = test4 - rSrad12;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_12 = dPP;

    if ( res2 < minres )
    {
        minres = res2;
        EVs_12 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_12 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_12 = dNN;
    }

    //std::cout << minres << std::endl;

    //Corner 22

    d1 = sqrt ( L2_22 ( 0 ) );
    d2 = sqrt ( L2_22 ( 1 ) );
    nd1 = -d1;
    nd2 = -d2;

    if ( d1 > 1.0 )
    {
        d1 = .98;
    }
    if ( d2 > 1.0 )
    {
        d2 = .98;
    }

    dPP = Matrix ( 2,2,d1,0.0,0.0,d2 );
    dPN = Matrix ( 2,2,d1,0.0,0.0,nd2 );
    dNP = Matrix ( 2,2,nd1,0.0,0.0,d2 );
    dNN = Matrix ( 2,2,nd1,0.0,0.0,nd2 );

    test1 = V1_22 * dPP * V2_22;
    test2 = V1_22 * dPN * V2_22;
    test3 = V1_22 * dNP * V2_22;
    test4 = V1_22 * dNN * V2_22;

    temp1 = test1 - rSrad22;
    temp2 = test2 - rSrad22;
    temp3 = test3 - rSrad22;
    temp4 = test4 - rSrad22;

    res1 = fabs ( temp1 ( 0,0 ) ) + fabs ( temp1 ( 0,1 ) ) + fabs ( temp1 ( 1,0 ) ) + fabs ( temp1 ( 1,1 ) );
    res2 = fabs ( temp2 ( 0,0 ) ) + fabs ( temp2 ( 0,1 ) ) + fabs ( temp2 ( 1,0 ) ) + fabs ( temp2 ( 1,1 ) );
    res3 = fabs ( temp3 ( 0,0 ) ) + fabs ( temp3 ( 0,1 ) ) + fabs ( temp3 ( 1,0 ) ) + fabs ( temp3 ( 1,1 ) );
    res4 = fabs ( temp4 ( 0,0 ) ) + fabs ( temp4 ( 0,1 ) ) + fabs ( temp4 ( 1,0 ) ) + fabs ( temp4 ( 1,1 ) );

    minres = res1;
    EVs_22 = dPP;



    if ( res2 < minres )
    {
        minres = res2;
        EVs_22 = dPN;
    }

    if ( res3 < minres )
    {
        minres = res3;
        EVs_22 = dNP;
    }

    if ( res4 < minres )
    {
        minres = res4;
        EVs_22 = dNN;
    }

//	EVs_11.print();
//	EVs_21.print();
//	EVs_12.print();
//	EVs_22.print();

    //std::cout << minres << std::endl;

//	std::cout << "****BEGIN QUAD****" << std::endl;
//
//	std::cout << "Corner 0,0: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_11.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_11.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_11.print();
//
//	std::cout << "Corner 1,0: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_21.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_21.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_21.print();
//
//	std::cout << "Corner 0,1: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_12.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_12.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_12.print();
//
//	std::cout << "Corner 1,1: " << std::endl;
//	std::cout << "Left:" << std::endl;
//	V1_22.getColumn(0).print();
//	std::cout << "Right: " << std::endl;
//	V2_22.getColumn(0).print();
//	std::cout << "Eigenvalues: " << std::endl;
//	EVs_22.print();
//	std::cout << "rSrad: " << std::endl;
//	rSrad22.print();
//
//	std::cout << "**** END QUAD ****" << std::endl;


    double Lam1_11, Lam2_11, Lam1_21, Lam2_21, Lam1_12, Lam2_12, Lam1_22, Lam2_22;
    Vector e1_11, e2_11, e1_21, e2_21, e1_12, e2_12, e1_22, e2_22;

    Lam1_11 = EVs_11 ( 0,0 );
    Lam2_11 = EVs_11 ( 1,1 );
    Lam1_21 = EVs_21 ( 0,0 );
    Lam2_21 = EVs_21 ( 1,1 );
    Lam1_12 = EVs_12 ( 0,0 );
    Lam2_12 = EVs_12 ( 1,1 );
    Lam1_22 = EVs_22 ( 0,0 );
    Lam2_22 = EVs_22 ( 1,1 );

    e1_11 = V1_11;
    e2_11 = V2_11;
    e1_21 = V1_21;
    e2_21 = V2_21;
    e1_12 = V1_12;
    e2_12 = V2_12;
    e1_22 = V1_22;
    e2_22 = V2_22;

    // Corner 1
//	if (L11(0) > L11(1))
//	{
//		Lam1_11 = L11(0);
//		Lam2_11 = L11(1);
//
//		e1_11 = V11.getColumn(0);
//		e2_11 = V11.getColumn(1);
//	}
//	else
//	{
//		Lam1_11 = L11(1);
//		Lam2_11 = L11(0);
//
//		e1_11 = V11.getColumn(1);
//		e2_11 = V11.getColumn(0);
//	}
//
//	// Corner 2
//	if (L21(0) > L21(1))
//	{
//		Lam1_21 = L21(0);
//		Lam2_21 = L21(1);
//
//		e1_21 = V21.getColumn(0);
//		e2_21 = V21.getColumn(1);
//	}
//	else
//	{
//		Lam1_21 = L21(1);
//		Lam2_21 = L21(0);
//
//		e1_21 = V21.getColumn(1);
//		e2_21 = V21.getColumn(0);
//	}
//
//	//Corner 3
//	if (L12(0) > L12(1))
//	{
//		Lam1_12 = L12(0);
//		Lam2_12 = L12(1);
//
//		e1_12 = V12.getColumn(0);
//		e2_12 = V12.getColumn(1);
//	}
//	else
//	{
//		Lam1_12 = L12(1);
//		Lam2_12 = L12(0);
//
//		e1_12 = V12.getColumn(1);
//		e2_12 = V12.getColumn(0);
//	}
//
//	//Corner 4
//	if (L22(0) > L22(1))
//	{
//		Lam1_22 = L22(0);
//		Lam2_22 = L22(1);
//
//		e1_22 = V22.getColumn(0);
//		e2_22 = V22.getColumn(1);
//	}
//	else
//	{
//		Lam1_22 = L22(1);
//		Lam2_22 = L22(0);
//
//		e1_22 = V22.getColumn(1);
//		e2_22 = V22.getColumn(0);
//	}

//	cout << endl << "Eigenvectors: " << endl;
//	cout << "(" << e1_11(0) << ", " << e1_11(1) << "), (" << e2_11(0) << ", " << e2_11(1) << ")" << endl;
//	cout << "(" << e1_21(0) << ", " << e1_21(1) << "), (" << e2_21(0) << ", " << e2_21(1) << ")" << endl;
//	cout << "(" << e1_12(0) << ", " << e1_12(1) << "), (" << e2_12(0) << ", " << e2_12(1) << ")" << endl;
//	cout << "(" << e1_22(0) << ", " << e1_22(1) << "), (" << e2_22(0) << ", " << e2_22(1) << ")" << endl;

    double testlambda[4] = { Lam1_11, 0, 0, Lam2_11 };

    Matrix testL = Matrix ( 2,2,testlambda,true );
    Matrix testinv;

    Matrix testsrad = V11 * testL * V11.inverse ( testinv );

//	cout << endl << endl;
//	cout << testsrad(0,0) << ", " << testsrad(0,1) << ", " << testsrad(1,0) << ", " << testsrad(1,1) << endl;
//	cout << endl << endl;

//	cout << endl << "Eigenvalues: " << endl;
//	cout << Lam1_11 << ", " << Lam2_11 << endl;
//	cout << Lam1_21 << ", " << Lam2_21 << endl;
//	cout << Lam1_12 << ", " << Lam2_12 << endl;
//	cout << Lam1_22 << ", " << Lam2_22 << endl;

    int num_iter = 10;
    double curru = 0;
    double currv = 0;
    double du = u / num_iter;
    double dv = v / num_iter;

    double logLam1_11 = log ( 1 - Lam1_11 );
    double logLam2_11 = log ( 1 - Lam2_11 );
    double logLam1_21 = log ( 1 - Lam1_21 );
    double logLam2_21 = log ( 1 - Lam2_21 );
    double logLam1_12 = log ( 1 - Lam1_12 );
    double logLam2_12 = log ( 1 - Lam2_12 );
    double logLam1_22 = log ( 1 - Lam1_22 );
    double logLam2_22 = log ( 1 - Lam2_22 );

    double logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    double logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    double avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    double avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    double Lam1 = 1 - exp ( logAvg1 );
    double Lam2 = 1 - exp ( logAvg2 );

    //Lam1 = logAvg1;
    //Lam2 = logAvg2;

    // Uses a form of bilinear interpolation to get the eigenvectors: This should probably be a Frechet mean of the thetas
    // on the unit circle in the future

    /*double theta1_11 = atan2(e1_11(1), e1_11(0));
    double theta2_11 = atan2(e2_11(1), e2_11(0));
    double theta1_21 = atan2(e1_21(1), e1_21(0));
    double theta2_21 = atan2(e2_21(1), e2_21(0));
    double theta1_12 = atan2(e1_12(1), e1_12(0));
    double theta2_12 = atan2(e2_12(1), e2_12(0));
    double theta1_22 = atan2(e1_22(1), e1_22(0));
    double theta2_22 = atan2(e2_22(1), e2_22(0));*/

    double PI = 3.14159265;

    double avgx1, avgy1, avgx2, avgy2, avgtheta1, avgtheta2;//, at1_11_21, at2_11_21, at1_12_22, at2_12_22;

//	std::cout << "eigs" << std::endl;
//	e1_11.print();
//	e1_21.print();
//	e1_12.print();
//	e1_22.print();

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    // Go from 11 to 21 (u direction)
    //double thetadist1_11_21 = theta1_21 - theta1_11;
    //if (thetadist1_11_21 > PI)
    //{
    //	theta1_21 -= 2*PI;
    //	thetadist1_11_21 = theta1_21 - theta1_11;
    //}
    //else if (thetadist1_11_21 < -PI)
    //{
    //	theta1_21 += 2*PI;
    //	thetadist1_11_21 = theta1_21 - theta1_11;
    //}
    //at1_11_21 = theta1_11 + u*thetadist1_11_21;

    //double thetadist2_11_21 = theta2_21 - theta2_11;
    //if (thetadist2_11_21 > PI)
    //{
    //	theta2_21 -= 2*PI;
    //	thetadist2_11_21 = theta2_21 - theta2_11;
    //}
    //else if (thetadist2_11_21 < -PI)
    //{
    //	theta2_21 += 2*PI;
    //	thetadist2_11_21 = theta2_21 - theta2_11;
    //}
    //at2_11_21 = theta2_11 + u*thetadist2_11_21;

    //// Go from 12 to 22 (u direction)
    //double thetadist1_12_22 = theta1_22 - theta1_12;
    //if (thetadist1_12_22 > PI)
    //{
    //	theta1_22 -= 2*PI;
    //	thetadist1_12_22 = theta1_22 - theta1_12;
    //}
    //else if (thetadist1_12_22 < -PI)
    //{
    //	theta1_22 += 2*PI;
    //	thetadist1_12_22 = theta1_22 - theta1_12;
    //}
    //at1_12_22 = theta1_12 + u*thetadist1_12_22;

    //double thetadist2_12_22 = theta2_22 - theta2_12;
    //if (thetadist2_12_22 > PI)
    //{
    //	theta2_22 -= 2*PI;
    //	thetadist2_12_22 = theta2_22 - theta2_12;
    //}
    //else if (thetadist2_12_22 < -PI)
    //{
    //	theta2_22 += 2*PI;
    //	thetadist2_12_22 = theta2_22 - theta2_12;
    //}
    //at2_12_22 = theta2_12 + u*thetadist2_12_22;

    //// Now go v direction
    //double thetadist1_v = at1_12_22 - at1_11_21;
    //if (thetadist1_v > PI)
    //{
    //	at1_12_22 -= 2*PI;
    //	thetadist1_v = at1_12_22 - at1_11_21;
    //}
    //else if (thetadist1_v < -PI)
    //{
    //	at1_12_22 += 2*PI;
    //	thetadist1_v = at1_12_22 - at1_11_21;
    //}
    //avgtheta1 = at1_11_21 + v*thetadist1_v;

    //double thetadist2_v = at2_12_22 - at2_11_21;
    //if (thetadist2_v > PI)
    //{
    //	at2_12_22 -= 2*PI;
    //	thetadist2_v = at2_12_22 - at2_11_21;
    //}
    //else if (thetadist2_v < -PI)
    //{
    //	at2_12_22 += 2*PI;
    //	thetadist2_v = at2_12_22 - at2_11_21;
    //}
    //avgtheta2 = at2_11_21 + v*thetadist2_v;

    double neweigenv[4] = { cos ( avgtheta1 ), sin ( avgtheta1 ), cos ( avgtheta2 ), sin ( avgtheta2 ) };
    Matrix newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    Matrix newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    Matrix NewV = Matrix ( 2,2,neweigenv,true );

    double newlambda[4] = { Lam1, 0, 0, Lam2 };

    Matrix NewL = Matrix ( 2,2,newlambda,true );

//	NewL.print();
    /*cout << V11(0,0) << ", " << V11(1,0)  << ", " << V11(0,1) << ", " << V11(1,1) << endl;
    cout << V21(0,0) << ", " << V21(1,0)  << ", " << V21(0,1) << ", " << V21(1,1) << endl;
    cout << V12(0,0) << ", " << V12(1,0)  << ", " << V12(0,1) << ", " << V12(1,1) << endl;
    cout << V22(0,0) << ", " << V22(1,0)  << ", " << V22(0,1) << ", " << V22(1,1) << endl;
    cout << NewV(0,0) << ", " << NewV(1,0) << ", " << NewV(0,1) << ", " << NewV(1,1) << endl;*/

    Matrix NewVi;
    NewV.inverse ( NewVi );

    //Matrix NewrSrad = NewV * NewL * NewVi;

    Matrix NewrSrad = newLeft * NewL * newRight;
//	NewrSrad.print();

    if ( side == 0 )
    {
        U0_11 = quadAtom11->getU0();
    }
    else
    {
        U0_11 = quadAtom11->getU1();
    }

    //Vector3D du11t_norm = du11t / du11t.norm();
    //Vector3D dv11t_norm = dv11t / dv11t.norm();



    double Pdata[6] = { du11t.getX(), dv11t.getX(), du11t.getY(), dv11t.getY(), du11t.getZ(), dv11t.getZ() };
    Matrix P = Matrix ( 2,3,Pdata,true );

//	du11.print();
//	dv11.print();
//	P.print();

    double Udata[3] = { U0_11.getX(), U0_11.getY(), U0_11.getZ() };
    Matrix U = Matrix ( 1,3,Udata,true );

    //double ru11 = du11 * U0_11;
    //double rv11 = dv11 * U0_11;

    //double ru21 = du21 * U0_21;
    //double rv21 = dv21 * U0_21;

    //double ru12 = du12 * U0_12;
    //double rv12 = dv12 * U0_12;

    //double ru22 = du22 * U0_22;
    //double rv22 = dv22 * U0_22;

//	cout << endl << "Ru & Rv: " << endl;
//	cout << ru11 << ", " << rv11 << endl;
//	cout << ru21 << ", " << rv21 << endl;
//	cout << ru12 << ", " << rv12 << endl;
//	cout << ru22 << ", " << rv22 << endl;

    double dSdu11data[6] = { dS0du11.getX(), dS0dv11.getX(), dS0du11.getY(), dS0dv11.getY(), dS0du11.getZ(), dS0dv11.getZ() };
    Matrix dSdu11 = Matrix ( 2,3,dSdu11data,true );

    double rdata[2] = { ru11, rv11 };
    Matrix rmat = Matrix ( 2,1,rdata,true );

    double Idata[9] = {1,0,0,0,1,0,0,0,1};
    Matrix I = Matrix ( 3,3,Idata,true );

    Matrix Q = P * ( ( U.t() * U ) - I );

    Matrix QQi;

    Matrix QQ = Q * Q.t();
    QQ.inverse ( QQi );

    Matrix UtU, UtUp, UtUi;

    UtU = U.t() * U;
	Q = P * ( UtU - I );
	UtUp = I + UtU;
	UtUp.inverse(UtUi);

	Matrix UtUpu = I + du*UtU;
	Matrix UtUpv = I + dv*UtU;

	Matrix UtUiu, UtUiv;
	UtUpu.inverse(UtUiu);
	UtUpv.inverse(UtUiv);


    Matrix R = -1 * P * U.t();

    Matrix dS = ( rSrad11.t() * Q ) + ( rmat * U );

    Matrix testrSrad11 = ( dSdu11 - ( rmat*U ) ) * ( Q.t() *QQi );
//	std::cout << "rs11: " << std::endl;
//	rSrad11.print();
//	std::cout << "trs11: " << std::endl;
//	testrSrad11.t().print();

    //Q.print();
    //Q11.print();

    Vector dSdu = dS.getRow ( 0 );
    Vector dSdv = dS.getRow ( 1 );

    //Vector3D dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
    //Vector3D dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );

    Vector Suu = ((rSrad11.t().getColumn(0) * Q.getRow(0)) -
    								2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
	Vector Suv = ((rSrad11.t().getColumn(1) * Q.getRow(1)) -
							2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

//		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
//						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
//		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
//						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

	Vector3D dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
	Vector3D dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

    Vector3D S0;

    if ( side == 0 )
    {
        S0 = U0_11 * quadAtom11->getR0();
    }
    else
    {
        S0 = U0_11 * quadAtom11->getR1();
    }


    Vector3D NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru = curru + du;
    currv = currv + dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    Vector3D newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );


    Matrix temprSrad = 0.25* ( rSrad11 + rSrad12 + rSrad21 + rSrad22 );

    //U0_11.print();
    for ( int i = 1; i < num_iter; i++ )
    {

//		std::cout << i << std::endl;
        // To find rSrad, need to calculate lambdas and eigenvectors
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        Q = P * ( ( U.t() * U ) - I );

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();

        UtU = U.t() * U;
		Q = P * ( UtU - I );
		UtUp = I + UtU;
		UtUp.inverse(UtUi);

		Matrix UtUpu = I + UtU;
		Matrix UtUpv = I + UtU;

		Matrix UtUiu, UtUiv;
		UtUpu.inverse(UtUiu);
		UtUpv.inverse(UtUiv);

		//Matrix Q2 = NewrSrad.t() * P



//        dS = ( NewrSrad.t() * Q ) + ( R * U );
//
//        dSdu = dS.getRow ( 0 );
//        dSdv = dS.getRow ( 1 );
//
//        dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
//        dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );
//
//        Vector3D TempSpoke = r * NewSpoke;
//
//        NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );
        //NewSpoke.print();

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
		Matrix S = Matrix(1,3,Spokedata,true);


		Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
								2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
		Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
								2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

//		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
//						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
//		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
//						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

		dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
		dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

		Vector3D TempSpoke = r*NewSpoke;
		NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);

        curru = curru + du;
        currv = currv + dv;

    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );

    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );
    Vector3D newbound = Vector3D ( Bxn ( 0,0 ), Byn ( 0,0 ), Bzn ( 0,0 ) );

    //cout << "{" << quadAtom11->getX().getX() + (quadAtom11->getR0() * quadAtom11->getU0().getX()) << ", " << quadAtom11->getX().getY() + (quadAtom11->getR0() * quadAtom11->getU0().getY()) << ", " << quadAtom11->getX().getZ() + (quadAtom11->getR0() * quadAtom11->getU0().getZ()) << "}" << endl;
    //cout << "{" << quadAtom21->getX().getX() + (quadAtom21->getR0() * quadAtom21->getU0().getX()) << ", " << quadAtom21->getX().getY() + (quadAtom21->getR0() * quadAtom21->getU0().getY()) << ", " << quadAtom21->getX().getZ() + (quadAtom21->getR0() * quadAtom21->getU0().getZ()) << "}" << endl;
    //cout << "{" << quadAtom12->getX().getX() + (quadAtom12->getR0() * quadAtom12->getU0().getX()) << ", " << quadAtom12->getX().getY() + (quadAtom12->getR0() * quadAtom12->getU0().getY()) << ", " << quadAtom12->getX().getZ() + (quadAtom12->getR0() * quadAtom12->getU0().getZ()) << "}" << endl;
    //cout << "{" << quadAtom22->getX().getX() + (quadAtom22->getR0() * quadAtom22->getU0().getX()) << ", " << quadAtom22->getX().getY() + (quadAtom22->getR0() * quadAtom22->getU0().getY()) << ", " << quadAtom22->getX().getZ() + (quadAtom22->getR0() * quadAtom22->getU0().getZ()) << "}" << endl;
    //cout << "{" << newpos.getX() + NewSpoke.getX() << ", " << newpos.getY() + NewSpoke.getY() << ", " << newpos.getZ() + NewSpoke.getZ() << endl;

    //cout << quadAtom11->getR0() << endl;
    //cout << quadAtom21->getR0() << endl;
    //cout << quadAtom12->getR0() << endl;
    //cout << quadAtom22->getR0() << endl;
    //cout << NewSpoke.norm() << endl;

    //cout << quadAtom11->getX().getX() << ", " << quadAtom11->getX().getY() << ", " << quadAtom11->getX().getZ() << endl;
    //cout << quadAtom21->getX().getX() << ", " << quadAtom21->getX().getY() << ", " << quadAtom21->getX().getZ() << endl;
    //cout << quadAtom12->getX().getX() << ", " << quadAtom12->getX().getY() << ", " << quadAtom12->getX().getZ() << endl;
    //cout << quadAtom22->getX().getX() << ", " << quadAtom22->getX().getY() << ", " << quadAtom22->getX().getZ() << endl;

    double NewR = NewSpoke.normalize();
    Vector3D htestD = newbound - newpos;
    double htestR = htestD.normalize();


    //cout << newpos.getX() << ", " << newpos.getY() << ", " << newpos.getZ() << endl;

    //cout << newpos.getX() << ", " << NewR << ", " << NewSpoke.getX() << endl;
    M3DSpoke hspoke = M3DSpoke ( newpos, htestD, htestR ); // This is for testing purposes
    //M3DSpoke* spoke11 = new M3DSpoke(newpos, htestD, htestR); // This is for testing purposes
    M3DSpoke spoke11 = M3DSpoke ( newpos, NewSpoke, NewR );

    // from corner 21

    curru = 1;
    currv = 0;
    du = ( u-1 ) /num_iter;
    dv = v / num_iter;

    //std::cout << "du: " << du << ", dv: " << dv << std::endl;

    logLam1_11 = log ( 1 - Lam1_11 );
    logLam2_11 = log ( 1 - Lam2_11 );
    logLam1_21 = log ( 1 - Lam1_21 );
    logLam2_21 = log ( 1 - Lam2_21 );
    logLam1_12 = log ( 1 - Lam1_12 );
    logLam2_12 = log ( 1 - Lam2_12 );
    logLam1_22 = log ( 1 - Lam1_22 );
    logLam2_22 = log ( 1 - Lam2_22 );

    logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    Lam1 = 1 - exp ( logAvg1 );
    Lam2 = 1 - exp ( logAvg2 );

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    neweigenv[0] = cos ( avgtheta1 );
    neweigenv[1] = sin ( avgtheta1 );
    neweigenv[2] = cos ( avgtheta2 );
    neweigenv[3] = sin ( avgtheta2 );

    newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    NewV = Matrix ( 2,2,neweigenv,true );

    newlambda[0] = Lam1;
    newlambda[1] = 0;
    newlambda[2] = 0;
    newlambda[3] = Lam2;

    NewL = Matrix ( 2,2,newlambda,true );

    NewV.inverse ( NewVi );

    NewrSrad = newLeft * NewL * newRight;


    Pdata[0] = du21t.getX();
    Pdata[1] = dv21t.getX();
    Pdata[2] = du21t.getY();
    Pdata[3] = dv21t.getY();
    Pdata[4] = du21t.getZ();
    Pdata[5] = dv21t.getZ();
    P = Matrix ( 2,3,Pdata,true );

    Udata[0] = U0_21.getX();
    Udata[1] = U0_21.getY();
    Udata[2] = U0_21.getZ();
    U = Matrix ( 1,3,Udata,true );

    double dSdu21data[6] = { dS0du21.getX(), dS0dv21.getX(), dS0du21.getY(), dS0dv21.getY(), dS0du21.getZ(), dS0dv21.getZ() };
    Matrix dSdu21 = Matrix ( 2,3,dSdu21data,true );

    rdata[0] = ru21;
    rdata[1] = rv21;
    rmat = Matrix ( 2,1,rdata,true );

    Q = P * ( ( U.t() * U ) - I );

    QQ = Q * Q.t();
    QQ.inverse ( QQi );

    R = -1 * P * U.t();

    dS = ( rSrad21.t() * Q ) + ( rmat * U );

    dSdu = dS.getRow ( 0 );
    dSdv = dS.getRow ( 1 );

    dSu = Vector3D ( dSdu ( 0 ),dSdu ( 1 ),dSdu ( 2 ) );
    dSv = Vector3D ( dSdv ( 0 ),dSdv ( 1 ),dSdv ( 2 ) );

    S0 = U0_21 * r21;

    NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru += du;
    currv += dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    for ( int i = 1; i < num_iter; i++ )
    {
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        double pdatam[6] = { npu.getX() - r*NewSpoke.getX(), npv.getX() - r*NewSpoke.getX(),
                		npu.getY() - r*NewSpoke.getY(), npv.getY() - r*NewSpoke.getY(),
                		npu.getZ() - r*NewSpoke.getZ(), npv.getZ() - r*NewSpoke.getZ()
                };
        Matrix Pm = Matrix(2,3,pdatam,true);
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        UtU = U.t() * U;
		Q = P * ( UtU - I );
		UtUp = I + UtU;
		UtUp.inverse(UtUi);

		Matrix UtUpu = I + UtU;
		Matrix UtUpv = I + UtU;

		Matrix UtUiu, UtUiv;
		UtUpu.inverse(UtUiu);
		UtUpv.inverse(UtUiv);

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();



        dS = ( NewrSrad.t() * Q ) + ( R * U );

        dSdu = dS.getRow ( 0 );
        dSdv = dS.getRow ( 1 );



        //Vector3D TempSpoke = r * NewSpoke;

        //NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
		Matrix S = Matrix(1,3,Spokedata,true);


		Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
								2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
		Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
								2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

//		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
//						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
//		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
//						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

		dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
		dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

        Vector3D TempSpoke = r*NewSpoke;
        NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);


        //NewSpoke.print();

        curru = curru + du;
        currv = currv + dv;
    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    NewR = NewSpoke.normalize();
    M3DSpoke spoke21 = M3DSpoke ( newpos, NewSpoke, NewR );

    // from corner 12

    curru = 0;
    currv = 1;
    du = u /num_iter;
    dv = ( v-1 ) / num_iter;

    //std::cout << "du: " << du << ", dv: " << dv << std::endl;

    logLam1_11 = log ( 1 - Lam1_11 );
    logLam2_11 = log ( 1 - Lam2_11 );
    logLam1_21 = log ( 1 - Lam1_21 );
    logLam2_21 = log ( 1 - Lam2_21 );
    logLam1_12 = log ( 1 - Lam1_12 );
    logLam2_12 = log ( 1 - Lam2_12 );
    logLam1_22 = log ( 1 - Lam1_22 );
    logLam2_22 = log ( 1 - Lam2_22 );

    logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    Lam1 = 1 - exp ( logAvg1 );
    Lam2 = 1 - exp ( logAvg2 );

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    neweigenv[0] = cos ( avgtheta1 );
    neweigenv[1] = sin ( avgtheta1 );
    neweigenv[2] = cos ( avgtheta2 );
    neweigenv[3] = sin ( avgtheta2 );

    newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    NewV = Matrix ( 2,2,neweigenv,true );

    newlambda[0] = Lam1;
    newlambda[3] = Lam2;

    NewL = Matrix ( 2,2,newlambda,true );

    NewV.inverse ( NewVi );

    NewrSrad = newLeft * NewL * newRight;

    Pdata[0] = du12t.getX();
    Pdata[1] = dv12t.getX();
    Pdata[2] = du12t.getY();
    Pdata[3] = dv12t.getY();
    Pdata[4] = du12t.getZ();
    Pdata[5] = dv12t.getZ();
    P = Matrix ( 2,3,Pdata,true );

    Udata[0] = U0_12.getX();
    Udata[1] = U0_12.getY();
    Udata[2] = U0_12.getZ();
    U = Matrix ( 1,3,Udata,true );

    double dSdu12data[6] = { dS0du12.getX(), dS0dv12.getX(), dS0du12.getY(), dS0dv12.getY(), dS0du12.getZ(), dS0dv12.getZ() };
    Matrix dSdu12 = Matrix ( 2,3,dSdu12data,true );

    rdata[0] = ru12;
    rdata[1] = rv12;
    rmat = Matrix ( 2,1,rdata,true );

    Q = P * ( ( U.t() * U ) - I );

    QQ = Q * Q.t();
    QQ.inverse ( QQi );

    R = -1 * P * U.t();

    dS = ( rSrad12.t() * Q ) + ( rmat * U );

    dSdu = dS.getRow ( 0 );
    dSdv = dS.getRow ( 1 );

    dSu = Vector3D ( dSdu ( 0 ),dSdu ( 1 ),dSdu ( 2 ) );
    dSv = Vector3D ( dSdv ( 0 ),dSdv ( 1 ),dSdv ( 2 ) );

    S0 = U0_12 * r12;

    NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru += du;
    currv += dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    for ( int i = 1; i < num_iter; i++ )
    {
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        Q = P * ( ( U.t() * U ) - I );

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();

        UtU = U.t() * U;
        		Q = P * ( UtU - I );
        		UtUp = I + UtU;
        		UtUp.inverse(UtUi);

        		Matrix UtUpu = I + UtU;
        		Matrix UtUpv = I + UtU;

        		Matrix UtUiu, UtUiv;
        		UtUpu.inverse(UtUiu);
        		UtUpv.inverse(UtUiv);

        dS = ( NewrSrad.t() * Q ) + ( R * U );

        dSdu = dS.getRow ( 0 );
        dSdv = dS.getRow ( 1 );

        dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
        dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );

        //Vector3D TempSpoke = r * NewSpoke;

        //NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );
        //NewSpoke.print();

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
        		Matrix S = Matrix(1,3,Spokedata,true);


        		Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
        								2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
        		Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
        								2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

        //		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
        //						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
        //		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
        //						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

        		dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
        		dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

                Vector3D TempSpoke = r*NewSpoke;
                NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);

        curru = curru + du;
        currv = currv + dv;
    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    NewR = NewSpoke.normalize();
    M3DSpoke spoke12 = M3DSpoke ( newpos, NewSpoke, NewR );

    // from corner 22

    curru = 1;
    currv = 1;
    du = ( u-1 ) /num_iter;
    dv = ( v-1 ) / num_iter;

    //std::cout << "du: " << du << ", dv: " << dv << std::endl;

    logLam1_11 = log ( 1 - Lam1_11 );
    logLam2_11 = log ( 1 - Lam2_11 );
    logLam1_21 = log ( 1 - Lam1_21 );
    logLam2_21 = log ( 1 - Lam2_21 );
    logLam1_12 = log ( 1 - Lam1_12 );
    logLam2_12 = log ( 1 - Lam2_12 );
    logLam1_22 = log ( 1 - Lam1_22 );
    logLam2_22 = log ( 1 - Lam2_22 );

    logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
    logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

    avg1 = ( 1-curru ) * ( 1-currv ) *Lam1_11 + ( curru ) * ( 1-currv ) *Lam1_21 + ( 1-curru ) * ( currv ) *Lam1_12 + ( curru ) * ( currv ) *Lam1_22;
    avg2 = ( 1-curru ) * ( 1-currv ) *Lam2_11 + ( curru ) * ( 1-currv ) *Lam2_21 + ( 1-curru ) * ( currv ) *Lam2_12 + ( curru ) * ( currv ) *Lam2_22;

    Lam1 = 1 - exp ( logAvg1 );
    Lam2 = 1 - exp ( logAvg2 );

    avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
    avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
    avgtheta1 = atan2 ( avgy1, avgx1 );

    avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
    avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
    avgtheta2 = atan2 ( avgy2, avgx2 );

    neweigenv[0] = cos ( avgtheta1 );
    neweigenv[1] = sin ( avgtheta1 );
    neweigenv[2] = cos ( avgtheta2 );
    neweigenv[3] = sin ( avgtheta2 );

    newLeft = Matrix ( 2,2,cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
    newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

    NewV = Matrix ( 2,2,neweigenv,true );

    newlambda[0] = Lam1;
    newlambda[3] = Lam2;

    NewL = Matrix ( 2,2,newlambda,true );

    NewV.inverse ( NewVi );

    NewrSrad = newLeft * NewL * newRight;


    Pdata[0] = du22t.getX();
    Pdata[1] = dv22t.getX();
    Pdata[2] = du22t.getY();
    Pdata[3] = dv22t.getY();
    Pdata[4] = du22t.getZ();
    Pdata[5] = dv22t.getZ();
    P = Matrix ( 2,3,Pdata,true );

    Udata[0] = U0_22.getX();
    Udata[1] = U0_22.getY();
    Udata[2] = U0_22.getZ();
    U = Matrix ( 1,3,Udata,true );

    double dSdu22data[6] = { dS0du22.getX(), dS0dv22.getX(), dS0du22.getY(), dS0dv22.getY(), dS0du22.getZ(), dS0dv22.getZ() };
    Matrix dSdu22 = Matrix ( 2,3,dSdu22data,true );

    rdata[0] = ru22;
    rdata[1] = rv22;
    rmat = Matrix ( 2,1,rdata,true );

    Q = P * ( ( U.t() * U ) - I );

    QQ = Q * Q.t();
    QQ.inverse ( QQi );

    R = -1 * P * U.t();

    dS = ( rSrad22.t() * Q ) + ( rmat * U );

    dSdu = dS.getRow ( 0 );
    dSdv = dS.getRow ( 1 );

    dSu = Vector3D ( dSdu ( 0 ),dSdu ( 1 ),dSdu ( 2 ) );
    dSv = Vector3D ( dSdv ( 0 ),dSdv ( 1 ),dSdv ( 2 ) );

    S0 = U0_22 * r22;

    NewSpoke = S0 + ( du * dSu ) + ( dv * dSv );

    curru += du;
    currv += dv;

    hu[0] = h1 ( curru );
    hu[1] = h2 ( curru );
    hu[2] = h3 ( curru );
    hu[3] = h4 ( curru );
    hv[0] = h1 ( currv );
    hv[1] = h2 ( currv );
    hv[2] = h3 ( currv );
    hv[3] = h4 ( currv );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    for ( int i = 1; i < num_iter; i++ )
    {
        logAvg1 = ( 1-curru ) * ( 1-currv ) *logLam1_11 + ( curru ) * ( 1-currv ) *logLam1_21 + ( 1-curru ) * ( currv ) *logLam1_12 + ( curru ) * ( currv ) *logLam1_22;
        logAvg2 = ( 1-curru ) * ( 1-currv ) *logLam2_11 + ( curru ) * ( 1-currv ) *logLam2_21 + ( 1-curru ) * ( currv ) *logLam2_12 + ( curru ) * ( currv ) *logLam2_22;

        //logAvg1 = (1-u)*(1-v)*logLam1_11 + (u)*(1-v)*logLam1_21 + (1-u)*(v)*logLam1_12 + (u)*(v)*logLam1_22;
        //logAvg2 = (1-u)*(1-v)*logLam2_11 + (u)*(1-v)*logLam2_21 + (1-u)*(v)*logLam2_12 + (u)*(v)*logLam2_22;
        Lam1 = 1 - exp ( logAvg1 );
        Lam2 = 1 - exp ( logAvg2 );

        avgx1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 0 ) + ( curru ) * ( 1-currv ) *e1_21 ( 0 ) + ( 1-curru ) * ( currv ) *e1_12 ( 0 ) + ( curru ) * ( currv ) *e1_22 ( 0 );
        avgy1 = ( 1-curru ) * ( 1-currv ) *e1_11 ( 1 ) + ( curru ) * ( 1-currv ) *e1_21 ( 1 ) + ( 1-curru ) * ( currv ) *e1_12 ( 1 ) + ( curru ) * ( currv ) *e1_22 ( 1 );
        avgtheta1 = atan2 ( avgy1, avgx1 );

        avgx2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 0 ) + ( curru ) * ( 1-currv ) *e2_21 ( 0 ) + ( 1-curru ) * ( currv ) *e2_12 ( 0 ) + ( curru ) * ( currv ) *e2_22 ( 0 );
        avgy2 = ( 1-curru ) * ( 1-currv ) *e2_11 ( 1 ) + ( curru ) * ( 1-currv ) *e2_21 ( 1 ) + ( 1-curru ) * ( currv ) *e2_12 ( 1 ) + ( curru ) * ( currv ) *e2_22 ( 1 );
        avgtheta2 = atan2 ( avgy2, avgx2 );

        neweigenv[0] = cos ( avgtheta1 );
        neweigenv[1] = sin ( avgtheta1 );
        neweigenv[2] = cos ( avgtheta2 );
        neweigenv[3] = sin ( avgtheta2 );
        NewV = Matrix ( 2,2,neweigenv,true );
        NewV.inverse ( NewVi );


        newLeft  = Matrix ( 2,2, cos ( avgtheta1 ), sin ( avgtheta1 ), sin ( avgtheta1 ), -cos ( avgtheta1 ) );
        newRight = Matrix ( 2,2, cos ( avgtheta2 ), sin ( avgtheta2 ), sin ( avgtheta2 ), -cos ( avgtheta2 ) );

        //newLeft.print();
        //newRight.print();
        newlambda[0] = Lam1;
        newlambda[1] = 0;
        newlambda[2] = 0;
        newlambda[3] = Lam2;
        NewL = Matrix ( 2,2,newlambda,true );
        //NewL.print();

        //NewrSrad = NewV * NewL * NewVi;
        NewrSrad = newLeft * NewL * newRight;
        //NewrSrad.print();

        double hu[4] = { h1 ( curru ), h2 ( curru ), h3 ( curru ), h4 ( curru ) };
        double hv[4] = { h1 ( currv ), h2 ( currv ), h3 ( currv ), h4 ( currv ) };
        Matrix humat = Matrix ( 1,4,hu,true );
        Matrix hvmat = Matrix ( 4,1,hv,true );

        double hup[4] = { h1p ( curru ), h2p ( curru ), h3p ( curru ), h4p ( curru ) };
        double hvp[4] = { h1p ( currv ), h2p ( currv ), h3p ( currv ), h4p ( currv ) };
        Matrix hupmat = Matrix ( 1,4,hup,true );
        Matrix hvpmat = Matrix ( 4,1,hvp,true );

        Matrix rup = hupmat * r_hermite_mat * hvmat;
        Matrix rvp = humat * r_hermite_mat * hvpmat;

        Matrix rp = humat * r_hermite_mat * hvmat;

        double rupd = rup ( 0,0 );
        double rvpd = rvp ( 0,0 );
        double rpd = rp ( 0,0 );

        Matrix pux = hupmat * hxmat * hvmat;
        Matrix puy = hupmat * hymat * hvmat;
        Matrix puz = hupmat * hzmat * hvmat;

        Matrix pvx = humat * hxmat * hvpmat;
        Matrix pvy = humat * hymat * hvpmat;
        Matrix pvz = humat * hzmat * hvpmat;

        Vector3D npu = Vector3D ( pux ( 0,0 ), puy ( 0,0 ), puz ( 0,0 ) );
        //npu.normalize();
        Vector3D npv = Vector3D ( pvx ( 0,0 ), pvy ( 0,0 ), pvz ( 0,0 ) );
        //npv.normalize();
        Vector3D nnorm = npu.cross ( npv );
        //nnorm.normalize();
        Vector3D npup = nnorm.cross ( npu );
        //npup.normalize();

        double pdata[6] = { npu.getX(), npv.getX(), npu.getY(), npv.getY(), npu.getZ(), npv.getZ() };
        P = Matrix ( 2,3,pdata,true );
        double r = NewSpoke.normalize();
        //std::cout << "r: " << r << std::endl;

        //r = rpd;
        //std::cout << "r2: " << rpd << std::endl;

        double udata[3] = { NewSpoke.getX(), NewSpoke.getY(), NewSpoke.getZ() };
        U = Matrix ( 1,3,udata,true );

        Q = P * ( ( U.t() * U ) - I );

        R = -1 * P * U.t();
        //R.print();
        R ( 0,0 ) = rupd;
        R ( 1,0 ) = rvpd;
        //R.print();

        UtU = U.t() * U;
        		Q = P * ( UtU - I );
        		UtUp = I + UtU;
        		UtUp.inverse(UtUi);

        		Matrix UtUpu = I + UtU;
        		Matrix UtUpv = I + UtU;

        		Matrix UtUiu, UtUiv;
        		UtUpu.inverse(UtUiu);
        		UtUpv.inverse(UtUiv);

        dS = ( NewrSrad.t() * Q ) + ( R * U );

        dSdu = dS.getRow ( 0 );
        dSdv = dS.getRow ( 1 );

        dSu = Vector3D ( dSdu ( 0 ), dSdu ( 1 ), dSdu ( 2 ) );
        dSv = Vector3D ( dSdv ( 0 ), dSdv ( 1 ), dSdv ( 2 ) );

        //Vector3D TempSpoke = r * NewSpoke;

        //NewSpoke = TempSpoke + ( du * dSu ) + ( dv * dSv );
        //NewSpoke.print();

        double Spokedata[3] = {r*NewSpoke.getX(), r*NewSpoke.getY(), r*NewSpoke.getZ()};
        		Matrix S = Matrix(1,3,Spokedata,true);


        		Vector Suu = ((NewrSrad.t().getColumn(0) * Q.getRow(0)) -
        								2*(P.getRow(0).dotProduct(U.getRow(0)) * U)) * UtUiu;
        		Vector Suv = ((NewrSrad.t().getColumn(1) * Q.getRow(1)) -
        								2*(P.getRow(1).dotProduct(U.getRow(0)) * U)) * UtUiv;

        //		dSdu = (NewrSrad.t().getColumn(0)*P.getRow(0)*Q) - (P.getRow(0) * UtU) -
        //						((P.getRow(0) - S.getRow(0))*UtU) - (Suu*UtU);
        //		dSdv = (NewrSrad.t().getColumn(1)*P.getRow(1)*Q) - (P.getRow(1) * UtU) -
        //						((P.getRow(1) - S.getRow(0))*UtU) - (Suv*UtU);

        		dSu = Vector3D ( Suu ( 0 ), Suu ( 1 ), Suu ( 2 ) );
        		dSv = Vector3D ( Suv ( 0 ), Suv ( 1 ), Suv ( 2 ) );

                Vector3D TempSpoke = r*NewSpoke;
                NewSpoke = TempSpoke + (du*dSu) + (dv*dSv);

        curru = curru + du;
        currv = currv + dv;
    }

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );

    NewR = NewSpoke.normalize();
    M3DSpoke spoke22 = M3DSpoke ( newpos, NewSpoke, NewR );

    // Mean of spokes:

    Vector3D s11 = spoke11.getS();
    Vector3D s21 = spoke21.getS();
    Vector3D s12 = spoke12.getS();
    Vector3D s22 = spoke22.getS();

    //spoke11->draw();
    //spoke12->draw();
    //spoke21->draw();
    //spoke22->draw();

    Vector3D meanspoke = ( ( 1-u ) * ( 1-v ) *s11 + ( u ) * ( 1-v ) *s21 + ( 1-u ) * ( v ) *s12 + ( u ) * ( v ) *s22 );
    NewR = meanspoke.normalize();

    M3DSpoke mean = M3DSpoke ( newpos, meanspoke, NewR );
    //mean->draw();

    //std::cout << "*****ENDING   INTERPOLATION*****" << std::endl;
    //return hspoke;
//         std::cout << "one" << std::endl;
    return mean;
    //return spoke11;

}

M3DSpoke M3DInterpolater::interpolateSpoke_hermite ( M3DFigure *figure, double u, double v, int side )
{
    int ubase = ( int ) floor ( u );
    int vbase = ( int ) floor ( v );

    u = u - ubase;
    v = v - vbase;

    M3DQuadFigure *tempFigure = dynamic_cast<M3DQuadFigure*> ( figure );

    if ( ubase == tempFigure->getRowCount() - 1 )
    {
        ubase = ubase - 1;
        u=1;
    }
    if ( vbase == tempFigure->getColumnCount() - 1 )
    {
        vbase = vbase - 1;
        v=1;
    }

    // Get four corner atoms

    M3DPrimitive *atom11 = tempFigure->getPrimitivePtr ( ubase,vbase );
    M3DQuadPrimitive* quadAtom11 = dynamic_cast<M3DQuadPrimitive*> ( atom11 );

    M3DPrimitive *atom21 = tempFigure->getPrimitivePtr ( ubase+1,vbase );
    M3DQuadPrimitive* quadAtom21 = dynamic_cast<M3DQuadPrimitive*> ( atom21 );

    M3DPrimitive *atom12 = tempFigure->getPrimitivePtr ( ubase,vbase+1 );
    M3DQuadPrimitive* quadAtom12 = dynamic_cast<M3DQuadPrimitive*> ( atom12 );

    M3DPrimitive *atom22 = tempFigure->getPrimitivePtr ( ubase+1,vbase+1 );
    M3DQuadPrimitive* quadAtom22 = dynamic_cast<M3DQuadPrimitive*> ( atom22 );

    Vector3D x11 = quadAtom11->getX();
    Vector3D x21 = quadAtom21->getX();
    Vector3D x12 = quadAtom12->getX();
    Vector3D x22 = quadAtom22->getX();

    double r11, r12, r21, r22, ru11, ru12, ru21, ru22, rv11, rv12, rv21, rv22;

    Vector3D U0_11, U0_12, U0_21, U0_22, B0_11, B0_12, B0_21, B0_22, S0_11, S0_12,S0_21,S0_22;

    if ( side == 0 )
    {
        U0_11 = quadAtom11->getU0();
        U0_21 = quadAtom21->getU0();
        U0_12 = quadAtom12->getU0();
        U0_22 = quadAtom22->getU0();

        B0_11 = quadAtom11->getX() + quadAtom11->getR0() * quadAtom11->getU0();
        B0_21 = quadAtom21->getX() + quadAtom21->getR0() * quadAtom21->getU0();
        B0_12 = quadAtom12->getX() + quadAtom12->getR0() * quadAtom12->getU0();
        B0_22 = quadAtom22->getX() + quadAtom22->getR0() * quadAtom22->getU0();

        S0_11 = U0_11 * quadAtom11->getR0();
        S0_21 = U0_21 * quadAtom21->getR0();
        S0_12 = U0_12 * quadAtom12->getR0();
        S0_22 = U0_22 * quadAtom22->getR0();

        r11 = quadAtom11->getR0();
        r21 = quadAtom21->getR0();
        r12 = quadAtom12->getR0();
        r22 = quadAtom22->getR0();

        ru11 = tempFigure->computeDR0du(ubase,vbase);
		ru12 = tempFigure->computeDR0du(ubase,vbase+1);
		ru21 = tempFigure->computeDR0du(ubase+1,vbase);
		ru22 = tempFigure->computeDR0du(ubase+1,vbase+1);

		rv11 = tempFigure->computeDR0dv(ubase,vbase);
		rv12 = tempFigure->computeDR0dv(ubase,vbase+1);
		rv21 = tempFigure->computeDR0dv(ubase+1,vbase);
		rv22 = tempFigure->computeDR0dv(ubase+1,vbase+1);

    }
    else
    {
        U0_11 = quadAtom11->getU1();
        U0_21 = quadAtom21->getU1();
        U0_12 = quadAtom12->getU1();
        U0_22 = quadAtom22->getU1();

        B0_11 = quadAtom11->getX() + quadAtom11->getR1() * quadAtom11->getU1();
        B0_21 = quadAtom21->getX() + quadAtom21->getR1() * quadAtom21->getU1();
        B0_12 = quadAtom12->getX() + quadAtom12->getR1() * quadAtom12->getU1();
        B0_22 = quadAtom22->getX() + quadAtom22->getR1() * quadAtom22->getU1();

        S0_11 = U0_11 * quadAtom11->getR1();
        S0_21 = U0_21 * quadAtom21->getR1();
        S0_12 = U0_12 * quadAtom12->getR1();
        S0_22 = U0_22 * quadAtom22->getR1();

        r11 = quadAtom11->getR1();
        r21 = quadAtom21->getR1();
        r12 = quadAtom12->getR1();
        r22 = quadAtom22->getR1();

        ru11 = tempFigure->computeDR1du(ubase,vbase);
		ru12 = tempFigure->computeDR1du(ubase,vbase+1);
		ru21 = tempFigure->computeDR1du(ubase+1,vbase);
		ru22 = tempFigure->computeDR1du(ubase+1,vbase+1);

		rv11 = tempFigure->computeDR1dv(ubase,vbase);
		rv12 = tempFigure->computeDR1dv(ubase,vbase+1);
		rv21 = tempFigure->computeDR1dv(ubase+1,vbase);
		rv22 = tempFigure->computeDR1dv(ubase+1,vbase+1);
    }

    Vector3D du11, du21, du12, du22, dv11, dv21, dv12, dv22;

    Vector3D dU0du11, dU0dv11, dU0du21, dU0dv21, dU0du12, dU0dv12, dU0du22, dU0dv22;
    Vector3D dB0du11, dB0dv11, dB0du21, dB0dv21, dB0du12, dB0dv12, dB0du22, dB0dv22;
    Vector3D dS0du11, dS0dv11, dS0du21, dS0dv21, dS0du12, dS0dv12, dS0du22, dS0dv22;

    dB0du11 = B0_21 - B0_11;
    dB0dv11 = B0_12 - B0_11;

    dB0du21 = B0_21 - B0_11;
    dB0dv21 = B0_22 - B0_21;

    dB0du12 = B0_22 - B0_12;
    dB0dv12 = B0_12 - B0_11;

    dB0du22 = B0_22 - B0_12;
    dB0dv22 = B0_22 - B0_21;

    du11 = x21 - x11;
    dv11 = x12 - x11;
    dU0du11 = U0_21 - U0_11;
    dU0dv11 = U0_12 - U0_11;
    dS0du11 = S0_21 - S0_11;
    dS0dv11 = S0_12 - S0_11;
    //ru11 = r21 - r11;
    //rv11 = r12 - r11;
//
//
    du21 = x21 - x11;
    dv21 = x22 - x21;
    dU0du21 = U0_21 - U0_11;
    dU0dv21 = U0_22 - U0_21;
    dS0du21 = S0_21 - S0_11;
    dS0dv21 = S0_22 - S0_21;
    //ru21 = r21 - r11;
    //rv21 = r22 - r21;
//
    du12 = x22 - x12;
    dv12 = x12 - x11;
    dU0du12 = U0_22 - U0_12;
    dU0dv12 = U0_12 - U0_11;
    dS0du12 = S0_22 - S0_12;
    dS0dv12 = S0_12 - S0_11;
    //ru12 = r22 - r12;
    //rv12 = r12 - r11;
//
    du22 = x22 - x12;
    dv22 = x22 - x21;
    dU0du22 = U0_22 - U0_12;
    dU0dv22 = U0_22 - U0_21;
    dS0du22 = S0_22 - S0_12;
    dS0dv22 = S0_22 - S0_21;
   // ru22 = r22 - r12;
  //  rv22 = r22 - r21;

    du11 = tempFigure->computeDXdu(ubase,vbase);
	du12 = tempFigure->computeDXdu(ubase,vbase+1);
	du21 = tempFigure->computeDXdu(ubase+1,vbase);
	du22 = tempFigure->computeDXdu(ubase+1,vbase+1);

	dv11 = tempFigure->computeDXdv(ubase,vbase);
	dv12 = tempFigure->computeDXdv(ubase,vbase+1);
	dv21 = tempFigure->computeDXdv(ubase+1,vbase);
	dv22 = tempFigure->computeDXdv(ubase+1,vbase+1);


    if ( ubase != 0 )
    {
        M3DQuadPrimitive* tempAtom = dynamic_cast<M3DQuadPrimitive*> ( tempFigure->getPrimitivePtr ( ubase-1, vbase ) );
    }


    // Get unit normals for the four corners and project derivatives on to the tangent plane
    Vector3D n11 = quadAtom11->getU0() - quadAtom11->getU1();
    n11.normalize();
    Vector3D n21 = quadAtom21->getU0() - quadAtom21->getU1();
    n21.normalize();
    Vector3D n12 = quadAtom12->getU0() - quadAtom12->getU1();
    n12.normalize();
    Vector3D n22 = quadAtom22->getU0() - quadAtom22->getU1();
    n22.normalize();

    Vector3D du11t = du11 - ( du11 * n11 ) * n11;
    Vector3D du21t = du21 - ( du21 * n21 ) * n21;
    Vector3D du12t = du12 - ( du12 * n12 ) * n12;
    Vector3D du22t = du22 - ( du22 * n22 ) * n22;

    Vector3D dv11t = dv11 - ( dv11 * n11 ) * n11;
    Vector3D dv21t = dv21 - ( dv21 * n21 ) * n21;
    Vector3D dv12t = dv12 - ( dv12 * n12 ) * n12;
    Vector3D dv22t = dv22 - ( dv22 * n22 ) * n22;

    // Build matrices for hermite interpolation of medial sheet
    double hx[16] = { x11.getX(), x21.getX(), du11t.getX(), du21t.getX(),
                      x12.getX(), x22.getX(), du12t.getX(), du22t.getX(), dv11t.getX(),
                      dv21t.getX(), 0, 0, dv12t.getX(), dv22t.getX(), 0, 0
                    };

    double hy[16] = { x11.getY(), x21.getY(), du11t.getY(), du21t.getY(),
                      x12.getY(), x22.getY(), du12t.getY(), du22t.getY(), dv11t.getY(),
                      dv21t.getY(), 0, 0, dv12t.getY(), dv22t.getY(), 0, 0
                    };

    double hz[16] = { x11.getZ(), x21.getZ(), du11t.getZ(), du21t.getZ(),
                      x12.getZ(), x22.getZ(),du12t.getZ(), du22t.getZ(), dv11t.getZ(),
                      dv21t.getZ(), 0, 0, dv12t.getZ(), dv22t.getZ(), 0, 0
                    };

    double r_hermite[16] = { r11, r21, ru11, ru21, r12, r22, ru12, ru22, rv11, rv21, 0, 0, rv12, rv22, 0, 0};

    Matrix r_hermite_mat = Matrix ( 4,4,r_hermite, true );

    Matrix hxmat = Matrix ( 4, 4, hx, true );
    Matrix hymat = Matrix ( 4, 4, hy, true );
    Matrix hzmat = Matrix ( 4, 4, hz, true );

    //U0_11.print();

    double hu[4] = { h1 ( u ), h2 ( u ), h3 ( u ), h4 ( u ) };
    double hv[4] = { h1 ( v ), h2 ( v ), h3 ( v ), h4 ( v ) };
    Matrix humat = Matrix ( 1, 4, hu, true );
    Matrix hvmat = Matrix ( 4, 1, hv, true );

    Matrix xn = humat * hxmat * hvmat;
    Matrix yn = humat * hymat * hvmat;
    Matrix zn = humat * hzmat * hvmat;

    double Bx[16] = { B0_11.getX(), B0_21.getX(), dB0du11.getX(), dB0du21.getX(),
                      B0_12.getX(), B0_22.getX(), dB0du12.getX(), dB0du22.getX(), dB0dv11.getX(),
                      dB0dv21.getX(), 0, 0, dB0dv12.getX(), dB0dv22.getX(), 0, 0
                    };

    double By[16] = { B0_11.getY(), B0_21.getY(), dB0du11.getY(), dB0du21.getY(),
                      B0_12.getY(), B0_22.getY(), dB0du12.getY(), dB0du22.getY(), dB0dv11.getY(),
                      dB0dv21.getY(), 0, 0, dB0dv12.getY(), dB0dv22.getY(), 0, 0
                    };

    double Bz[16] = { B0_11.getZ(), B0_21.getZ(), dB0du11.getZ(), dB0du21.getZ(),
                      B0_12.getZ(), B0_22.getZ(), dB0du12.getZ(), dB0du22.getZ(), dB0dv11.getZ(),
                      dB0dv21.getZ(), 0, 0, dB0dv12.getZ(), dB0dv22.getZ(), 0, 0
                    };

    Matrix bxmat = Matrix ( 4, 4, Bx, true );
    Matrix bymat = Matrix ( 4, 4, By, true );
    Matrix bzmat = Matrix ( 4, 4, Bz, true );

    Matrix Bxn = humat * bxmat * hvmat;
    Matrix Byn = humat * bymat * hvmat;
    Matrix Bzn = humat * bzmat * hvmat;

    hu[0] = h1 ( u );
    hu[1] = h2 ( u );
    hu[2] = h3 ( u );
    hu[3] = h4 ( u );
    hv[0] = h1 ( v );
    hv[1] = h2 ( v );
    hv[2] = h3 ( v );
    hv[3] = h4 ( v );

    humat = Matrix ( 1,4,hu,true );
    hvmat = Matrix ( 4,1,hv,true );

    xn = humat * hxmat * hvmat;
    yn = humat * hymat * hvmat;
    zn = humat * hzmat * hvmat;

    Vector3D newpos = Vector3D ( xn ( 0,0 ), yn ( 0,0 ), zn ( 0,0 ) );
    Vector3D newbound = Vector3D ( Bxn ( 0,0 ), Byn ( 0,0 ), Bzn ( 0,0 ) );


    Vector3D htestD = newbound - newpos;
    double htestR = htestD.normalize();

    M3DSpoke hspoke = M3DSpoke ( newpos, htestD, htestR ); // This is for testing purposes

    return hspoke;

}
void M3DInterpolater::interpolate2D(vnl_double_3x1& x1, vnl_double_3x1& x2, double d, vnl_double_3x1& dx1, vnl_double_3x1& dx2, vnl_double_3x1& x_interp)
{
    // distance from left-most point
//    double d = 0.5;

    // hermite coefficients
    double coef1 = h1(d);
    double coef2 = h2(d);
    double coef3 = h3(d);
    double coef4 = h4(d);

    // this is skeletal point pos
    x_interp = coef1 * x1 + coef2* x2 + coef3 * dx1 + coef4 * dx2;

}
void M3DInterpolater::interpolate2D(double r1, double r2, vnl_double_3x1& U1, vnl_double_3x1& U2, vnl_double_3x1& x1, vnl_double_3x1& x2, vnl_double_3x1& dx1, vnl_double_3x1& dx2, int level, double lamda, std::vector<M3DSpoke>& output)
{
    if(level == 0)
    {
        return;
    }
    double d = lamda/2;

    // 1. interpolate skeletal points for each step
    vnl_double_3x1 x_interp;
    interpolate2D(x1, x2, d, dx1, dx2, x_interp);

    // 2. interpolate U and r
    vnl_double_3x1 uMid;
    double rMid = computeMiddleS(U1, U2, r1*U1, r2*U2, lamda, uMid);
    output.push_back(M3DSpoke(x_interp(0,0), x_interp(1,0), x_interp(2,0), uMid(0,0), uMid(1,0), uMid(2,0), rMid));

    // 3. interpolate finer
    // left side of middle point
    interpolate2D(r1, rMid, U1, uMid, x1, x_interp, dx1, dx2, level-1, lamda/2, output);
    // right side of middle point
    interpolate2D(rMid, r2, uMid, U2, x_interp, x2, dx1, dx2, level-1, lamda/2, output);

    return;
}

void M3DInterpolater::printInorder(int level, std::vector<M3DSpoke>& input, std::vector<M3DSpoke>& output)
{
    // if(level == 0)
    // {
    //     return;
    // }

    if(input.size()==1)
    {
        output.push_back(input[0]);
        return;
    }
    int midIndex = input.size() / 2;
    std::vector<M3DSpoke> leftTree, rightTree;
    for(int i = 1; i < midIndex + 1; ++i)
    {
        leftTree.push_back(input[i]);
    }
    for(int i = midIndex+1; i < input.size(); ++i)
    {
        rightTree.push_back(input[i]);
    }
    // print left tree
    printInorder(level-1, leftTree, output);

    // print middle, for the current input, the first is always middle spoke
    output.push_back(input[0]);

    // print right tree
    printInorder(level-1, rightTree, output);
}
// main entry of interpolate 2D srep
void M3DInterpolater::interpolate2D(M3DFigure *figure, int level, int side, std::vector<M3DSpoke>& output)
{
    output.clear();
    M3DQuadFigure* fig = dynamic_cast<M3DQuadFigure*>(figure);
    int primitiveCount = figure->getPrimitiveCount();

    // add first crest spoke
    output.push_back(figure->getPrimitivePtr(0)->getSpoke(2));
    for(int i = 0; i < primitiveCount-1; ++i)
    {
        std::vector<M3DSpoke> interpolateResult;
        M3DPrimitive* currentPrimitive = figure->getPrimitivePtr(i);
        M3DPrimitive* nextPrimitive = figure->getPrimitivePtr(i+1);
        vnl_double_3x1 x1, x2;

        convertVector3DToVNL(currentPrimitive->getX(), x1);
        convertVector3DToVNL(nextPrimitive->getX(), x2);

        vnl_double_3x1 dx1, dx2;
        convertVector3DToVNL(fig->computeDX(i), dx1);
        convertVector3DToVNL(fig->computeDX(i+1), dx2);

        vnl_double_3x1 U1, U2;
        double r1, r2;

        r1 = getR(currentPrimitive, side);
        r2 = getR(nextPrimitive, side);

        getU(currentPrimitive, side, U1);
        getU(nextPrimitive, side, U2);
//        interpolate2D(r1, r2, U1, U2, x1, x2, dx1, dx2, level, 1, output);
        interpolate2D(r1, r2, U1, U2, x1, x2, dx1, dx2, level, 1, interpolateResult);
        output.push_back(M3DSpoke(x1(0,0), x1(1,0), x1(2,0), U1(0,0), U1(1,0), U1(2,0), r1));

        // adjust the order of spokes along the skeleton, like tranverse a tree in inorder
        printInorder(level, interpolateResult, output);

    }

    // add last normal spoke
    M3DPrimitive* lastPrimitive = figure->getPrimitivePtr(primitiveCount - 1);
    output.push_back(lastPrimitive->getSpoke(side));

    // add last crest spoke
    output.push_back(lastPrimitive->getSpoke(2));

    // interpolate crest spoke
}
void  M3DInterpolater::getSpokeFromInterpolation(M3DFigure* figure, int side, double x, M3DSpoke& output)
{
    // decide interpolate level
    int interpolateLevel = 0;
    double resolution = 1.0; // the finest step size in parameter space
    double fraction = x - int(x); // the fraction part of x

    // when the fraction of x is multiple of resolution, stop dividing
    // i.e. when factor is integer
//
    double factor = fraction / resolution;
    double test = 10.0;
    do
    {
        test = ((fraction/resolution) - int(fraction/resolution));
        resolution /= 2;
        interpolateLevel++;
        if((test < 1e-4 && test >= 0) || resolution < 1e-4 || interpolateLevel > 8) break;
    }while(true);

    factor = int(fraction / resolution);

    // interpolate all spokes
    std::vector<M3DSpoke> spokes;
    this->interpolate2D(figure, interpolateLevel, side, spokes);

    // select the spoke at position x
    int numInInterval = pow(2, interpolateLevel) - 1; // the number of spokes in each interval
    int index = int(x) + int(x) * numInInterval + factor + 1; // the first spoke in spokes is crest spoke, thus need +1
    output.setR(spokes[index].getR());
    output.setU(spokes[index].getU());
    output.setX(spokes[index].getX());
}