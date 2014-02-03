#ifndef CUBICHERMITESPLINE_H
#define CUBICHERMITESPLINE_H

#include <unordered_map>

#include "spline_library/spline.h"

class CubicHermiteSpline : public Spline
{
//constructors
public:
    CubicHermiteSpline(const std::vector<Vector3D> &points, const std::vector<Vector3D> &tangents, double alpha = 0.0);
protected:
    //you're only allowed to create one of these without point data if a subclass is providing the point data
    CubicHermiteSpline();

//methods
public:
    virtual Vector3D getPosition(double x) const;
    virtual InterpolatedPT getTangent(double x) const;
    virtual InterpolatedPTC getCurvature(double x) const;
    virtual InterpolatedPTCW getWiggle(double x) const;

    virtual double getT(int index) const;
    virtual double getMaxT(void) const;

    virtual const std::vector<Vector3D> &getPoints(void) const;

    virtual bool isLooping(void) const;

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

struct CubicHermiteSpline::InterpolationData
{
    //points
    Vector3D p0;
    Vector3D p1;

    //tangents
    Vector3D m0;
    Vector3D m1;

    //t values
    double t0;
    double t1;

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

    //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
    //intuitively it would just be the derivative of the position function and nothing else
    //if you know why please let me know
    return (
            d_basis00 * segment.p0 +
            d_basis10 * segment.m0 +

            d_basis11 * segment.m1 +
            d_basis01 * segment.p1
            ) * segment.tDistanceInverse;
}

inline Vector3D CubicHermiteSpline::computeCurvature(double t, const InterpolationData &segment) const
{
    //calculate the acceleration of the spline at t.
    //ie just compute the second derivatives of all the basis functions
    double d2_basis00 = 6 * (2 * t - 1);
    double d2_basis10 = 2 * (3 * t - 2);

    double d2_basis11 = 2 * (3 * t - 1);
    double d2_basis01 = -d2_basis00;

    //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
    //intuitively it would just be the 2nd derivative of the position function and nothing else
    //if you know why please let me know
    return (
            d2_basis00 * segment.p0 +
            d2_basis10 * segment.m0 +

            d2_basis11 * segment.m1 +
            d2_basis01 * segment.p1
            ) * (segment.tDistanceInverse * segment.tDistanceInverse);
}

inline Vector3D CubicHermiteSpline::computeWiggle(double t, const InterpolationData &segment) const
{
    //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
    //intuitively it would just be the 2nd derivative of the position function and nothing else
    //if you know why please let me know
    return (
                12 * (segment.p0 - segment.p1) + 6 * (segment.m0 + segment.m1)
            ) * (segment.tDistanceInverse * segment.tDistanceInverse * segment.tDistanceInverse);
}

#endif // CUBICHERMITESPLINE_H
