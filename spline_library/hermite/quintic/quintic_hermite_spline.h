#ifndef QUINTICHERMITESPLINE_H
#define QUINTICHERMITESPLINE_H

#include <unordered_map>

#include "spline_library/spline.h"

class QuinticHermiteSpline : public Spline
{
//constructors
public:
    QuinticHermiteSpline(const std::vector<Vector3D> &points,
                         const std::vector<Vector3D> &tangents,
                         const std::vector<Vector3D> &curvatures,
                         float alpha
                         );
protected:
    //you're only allowed to create one of these without point data if a subclass is providing the point data
    QuinticHermiteSpline();

//methods
public:
    virtual Vector3D getPosition(double x) const override;
    virtual InterpolatedPT getTangent(double x) const override;
    virtual InterpolatedPTC getCurvature(double x) const override;
    virtual InterpolatedPTCW getWiggle(double x) const override;

    virtual double getT(int index) const override;
    virtual double getMaxT(void) const override;

    virtual const std::vector<Vector3D> &getPoints(void) const override;

    virtual bool isLooping(void) const override;

protected:
	struct InterpolationData;

    inline Vector3D computePosition(double t, const InterpolationData &segment) const;
	inline Vector3D computeTangent(double t, const InterpolationData &segment) const;
    inline Vector3D computeCurvature(double t, const InterpolationData &segment) const;
    inline Vector3D computeWiggle(double t, const InterpolationData &segment) const;

	int getSegmentIndex(double x) const;

//data
protected:
	//a vector containing pre-computed datasets, one per segment
    //there will be lots of duplication of data here,
    //but precomputing this really speeds up the interpolation
	int numSegments;
    std::vector<InterpolationData> segmentData;

    double maxT;

    //original point data
    std::vector<Vector3D> points;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,double> indexToT;
};

struct QuinticHermiteSpline::InterpolationData
{
    //points
    Vector3D p0;
    Vector3D p1;

    //tangents
    Vector3D m0;
    Vector3D m1;

    //curvatures
    Vector3D c0;
    Vector3D c1;

    //t values
    double t0;
    double t1;

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
    
    return (
		basis00 * segment.p0 + 
		basis10 * segment.m0 + 
		basis20 * segment.c0 +
		basis21 * segment.c1 + 
		basis11 * segment.m1 +
        basis01 * segment.p1
            ) * segment.tDistanceInverse;
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

    return (
        basis00 * segment.p0 +
        basis10 * segment.m0 +
        basis20 * segment.c0 +
        basis21 * segment.c1 +
        basis11 * segment.m1 +
        basis01 * segment.p1
            ) * (segment.tDistanceInverse * segment.tDistanceInverse);
}

inline Vector3D QuinticHermiteSpline::computeWiggle(double t, const InterpolationData &segment) const
{
    //we're essentially computing the third derivative of the computePosition function with respect to t
    double basis00 = -60 * (6 * t * (t - 1) + 1);
    double basis10 = -12 * (t * (15 * t - 16) + 3);
    double basis20 = t * (36 - 30 * t) - 9;
    double basis21 = t * (30 * t - 24) + 3;
    double basis11 = -12 * (t * (15 * t - 14) + 2);
    double basis01 = -basis00;

    return (
        basis00 * segment.p0 +
        basis10 * segment.m0 +
        basis20 * segment.c0 +
        basis21 * segment.c1 +
        basis11 * segment.m1 +
        basis01 * segment.p1
            ) * (segment.tDistanceInverse * segment.tDistanceInverse * segment.tDistanceInverse);
}

#endif // QUINTICHERMITESPLINE_H
