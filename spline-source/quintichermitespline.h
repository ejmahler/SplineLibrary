#ifndef QUINTICHERMITESPLINE_H
#define QUINTICHERMITESPLINE_H

#include "spline.h"

class QuinticHermiteSpline : public Spline
{
protected: //methods
	struct InterpolationData;

    inline Vector3D computePosition(double t, const InterpolationData &segment) const;
    Vector3D computePosition2(double t, const InterpolationData &segment) const;
	inline Vector3D computeTangent(double t, const InterpolationData &segment) const;
	inline Vector3D computeCurvature(double t, const InterpolationData &segment) const;

	int getSegmentIndex(double x) const;

protected: //data
	//a vector containing pre-computed datasets, one per segment
    //there will be lots of duplication of data here,
    //but precomputing this really speeds up the interpolation
	int numSegments;
    std::vector<InterpolationData> segmentData;


};

struct QuinticHermiteSpline::InterpolationData {
    //t values
    double t0;
    double t1;

    //points
    Vector3D p0;
    Vector3D p1;

    //tangents
    Vector3D m0;
    Vector3D m1;

    //curvatures
    Vector3D c0;
    Vector3D c1;

    //reciprocal of distance in T between p0 and p1
    double tDistanceInverse;

    //padding to make sure we're aligned to 16 bytes
    double padding;
};

inline Vector3D QuinticHermiteSpline::computePosition(double t, const InterpolationData &segment) const
{
    //interpolate the position of the given segment at t

    double oneMinusT = 1 - t;

    //this is a logical extension of the cubic hermite spline's basis functions
    //that has one basis function for t0 position, one for t1 position
    //one for t0 tangent (1st derivative of position), and one for t1 tangent
    //this adds 2 more basis functions, one for t0 curvature (2nd derivative) and t1 curvature
    //see this paper for details http://www.rose-hulman.edu/~finn/CCLI/Notes/day09.pdf
    double basis00 = (oneMinusT * oneMinusT * oneMinusT * (t * (6 * t + 3) + 1));
    double basis10 = (t * oneMinusT * oneMinusT * oneMinusT * (3 * t + 1));
    double basis20 = 0.5 * oneMinusT * oneMinusT * oneMinusT * t * t;
    double basis21 = 0.5 * oneMinusT * oneMinusT * t * t * t;
    double basis11 = t * t * t * oneMinusT * (t * 3 - 4);
    double basis01 = t * t * t * (t * (6 * t - 15) + 10);

    return
        basis00 * segment.p0 +
        basis10 * segment.m0 +
        basis20 * segment.c0 +
        basis21 * segment.c1 +
        basis11 * segment.m1 +
        basis01 * segment.p1;
}

inline Vector3D QuinticHermiteSpline::computeTangent(double t, const InterpolationData &segment) const
{
	//interpolate the first derivative of the given segment at t

	double oneMinusT = 1 - t;

	//we're essentially computing the derivative of the computePosition function with respect to t
	//we can do this by computing the derivatives of each of its basis functions.
	//thankfully this can easily be done analytically since they're polynomials!
    double basis00 = -30 * oneMinusT * oneMinusT * t * t;
    double basis10 = -oneMinusT * oneMinusT * (3 * t - 1) * (5 * t + 1);
	double basis20 = -0.5 * oneMinusT * oneMinusT * t * (5 * t - 2);
	double basis21 = -0.5 * oneMinusT * t * t * (5 * t - 3);
	double basis11 = -t * t * (3 * t - 2) * (5 * t - 6);
    double basis01 = 30 * oneMinusT * oneMinusT * t * t;
    
    return 
		basis00 * segment.p0 + 
		basis10 * segment.m0 + 
		basis20 * segment.c0 +
		basis21 * segment.c1 + 
		basis11 * segment.m1 +
		basis01 * segment.p1;
}

inline Vector3D QuinticHermiteSpline::computeCurvature(double t, const InterpolationData &segment) const
{
	//interpolate the second derivative of the given segment at t

	double oneMinusT = 1 - t;

	//we're essentially computing the second derivative of the computePosition function with respect to t
	//we can do this by computing the second derivatives of each of its basis functions.
	//thankfully this can easily be done analytically since they're polynomials!
    double basis00 = 60 * oneMinusT * t * (2 * t - 1);
    double basis10 = 12 * oneMinusT * t * (5 * t - 3);
	double basis20 = t * (t * (-10 * t + 18) - 9) + 1;
	double basis21 = t * (t * (10 * t - 12) + 3);
	double basis11 = 12 * oneMinusT * t * (5 * t - 2);
    double basis01 = -60 * oneMinusT * t * (2 * t - 1);
    
    return 
		basis00 * segment.p0 + 
		basis10 * segment.m0 + 
		basis20 * segment.c0 +
		basis21 * segment.c1 + 
		basis11 * segment.m1 +
		basis01 * segment.p1;
}

#endif // QUINTICHERMITESPLINE_H
