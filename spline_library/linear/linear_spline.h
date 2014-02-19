#ifndef LINEAR_SPLINE_H
#define LINEAR_SPLINE_H

#include <vector>
#include <unordered_map>

#include "spline_library/spline.h"

class LinearSpline : public Spline
{
//constructors
public:
    LinearSpline(const std::vector<Vector3D> &points, double alpha = 0.0);
protected:
    //you're only allowed to create one of these without point data if a subclass is providing the point data
    LinearSpline();

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
    inline Vector3D computeTangent(const InterpolationData &segment) const;

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

struct LinearSpline::InterpolationData
{
    //points
    Vector3D p0;
    Vector3D p1;

    //t values
    double t0;
    double t1;

    //reciprocal of distance in T between p0 and p1
    double tDistanceInverse;

    //padding to make sure we're aligned to 16 bytes
    double padding;
};

inline Vector3D LinearSpline::computePosition(double t, const InterpolationData &segment) const
{
    double oneMinusT = 1 - t;

    return segment.p0 * oneMinusT + segment.p1 * t;
}

inline Vector3D LinearSpline::computeTangent(const InterpolationData &segment) const
{
    return (segment.p1 - segment.p0) * segment.tDistanceInverse;
}

#endif // LINEAR_SPLINE_H
