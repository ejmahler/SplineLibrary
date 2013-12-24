#ifndef CUBICHERMITESPLINE_H
#define CUBICHERMITESPLINE_H

#include <unordered_map>

#include "../spline.h"

class CubicHermiteSpline : public Spline
{
public:
    virtual Vector3D getPosition(double x) const;
    virtual InterpolatedPT getTangent(double x) const;
    virtual InterpolatedPTC getCurvature(double x) const;

    virtual double getT(int index) const;
    virtual double getMaxT(void) const;
    virtual int getNumSegments(void) const;

    virtual const std::vector<Vector3D> &getPoints(void) const;

    virtual bool isLooping(void) const;

protected: //methods
    struct InterpolationData;

    inline Vector3D computePosition(double t, const InterpolationData &segment) const;
    inline Vector3D computeTangent(double t, const InterpolationData &segment) const;
    inline Vector3D computeCurvature(double t, const InterpolationData &segment) const;

    int getSegmentIndex(double x) const;

protected: //data
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

struct CubicHermiteSpline::InterpolationData {
    //t values
    double t0;
    double t1;

    //points
    Vector3D p0;
    Vector3D p1;

    //tangents
    Vector3D m0;
    Vector3D m1;

    //reciprocal of distance in T between p0 and p1
    double tDistanceInverse;

    //padding to make sure we're aligned to 16 bytes
    double padding;
};

inline Vector3D CubicHermiteSpline::computePosition(double t, const InterpolationData &segment) const
{
    double oneMinusT = 1 - t;

    double basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
    double basis10 = t * oneMinusT * oneMinusT;

    double basis11 = t * t * -oneMinusT;
    double basis01 = t * t * (3 - 2*t);

    return
            basis00 * segment.p0 +
            basis10 * segment.m0 +

            basis11 * segment.m1 +
            basis01 * segment.p1;
}

inline Vector3D CubicHermiteSpline::computeTangent(double t, const InterpolationData &segment) const
{
    double oneMinusT = 1 - t;

    //calculate the velocity of the spline at t.
    //ie just compute the derivatives of all the basis functions
    double d_basis00 = -6 * t * oneMinusT;
    double d_basis10 = -(3*t - 1) * oneMinusT;

    double d_basis11 = t * (3 * t - 2);
    double d_basis01 = -d_basis00;

    return
            d_basis00 * segment.p0 +
            d_basis10 * segment.m0 +

            d_basis11 * segment.m1 +
            d_basis01 * segment.p1;
}

inline Vector3D CubicHermiteSpline::computeCurvature(double t, const InterpolationData &segment) const
{
    //calculate the acceleration of the spline at t.
    //ie just compute the second derivatives of all the basis functions
    double d2_basis00 = 6 * (2 * t - 1);
    double d2_basis10 = 2 * (3 * t - 2);

    double d2_basis11 = 2 * (3 * t - 1);
    double d2_basis01 = -d2_basis00;

    return
            d2_basis00 * segment.p0 +
            d2_basis10 * segment.m0 +

            d2_basis11 * segment.m1 +
            d2_basis01 * segment.p1;
}

#endif // CUBICHERMITESPLINE_H
