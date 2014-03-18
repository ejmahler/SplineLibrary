#ifndef B_SPLINE_H
#define B_SPLINE_H

#include "spline_library/spline.h"

#include <unordered_map>

class CubicBSpline : public Spline
{
//constructors
public:
    CubicBSpline(const std::vector<Vector3D> &points);
protected:
    //you're only allowed to create one of these without point data if a subclass is providing the point data
    CubicBSpline();

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
    inline Vector3D computeWiggle(const InterpolationData &segment) const;

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

struct CubicBSpline::InterpolationData {
    //t values
    double t1, t2;

    //points
    Vector3D p0, p1, p2, p3;
};

inline Vector3D CubicBSpline::computePosition(double t, const InterpolationData &segment) const
{
    return (
                segment.p0 * ((1 - t) * (1 - t) * (1 - t)) +
                segment.p1 * (t * t * 3 * (t - 2) + 4) +
                segment.p2 * (t * (t * (-3 * t + 3) + 3) + 1) +
                segment.p3 * (t * t * t)
            ) / 6;
}

inline Vector3D CubicBSpline::computeTangent(double t, const InterpolationData &segment) const
{
    return (
                segment.p0 * (-(1 - t) * (1 - t)) +
                segment.p1 * (t * (3 * t - 4)) +
                segment.p2 * ((3 * t + 1) * (1 - t)) +
                segment.p3 * (t * t)
            ) / 2;
}

inline Vector3D CubicBSpline::computeCurvature(double t, const InterpolationData &segment) const
{
    return (
                segment.p0 * (1 - t) +
                segment.p1 * (3 * t - 2) +
                segment.p2 * (1 - 3 * t) +
                segment.p3 * (t)
            );
}

inline Vector3D CubicBSpline::computeWiggle(const InterpolationData &segment) const
{
    return 3 * (segment.p1 - segment.p2) + (segment.p3 - segment.p0);
}

#endif // B_SPLINE_H
