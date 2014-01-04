#include "looping_cubic_hermite_spline.h"

#include <cmath>
#include <cassert>

LoopingCubicHermiteSpline::LoopingCubicHermiteSpline()
{
}

LoopingCubicHermiteSpline::LoopingCubicHermiteSpline(const std::vector<Vector3D> &points, const std::vector<Vector3D> &tangents, double alpha)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());

    this->points = points;

    std::unordered_map<int, double> indexToT_Raw;

    int size = points.size();
    numSegments = size;

    //we know points[0] will have a t value of 0
    indexToT_Raw[0] = 0;

    //compute the t values of the other points
    for(int i = 0; i < size + 1; i++)
    {
        double distance = (points.at(i%size) - points.at((i - 1)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing tem by maxT
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = numSegments * it->second / maxTRaw;
    }
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangents.at(i) * tDistance;
        segment.m1 = tangents.at(i + 1) * tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D LoopingCubicHermiteSpline::getPosition(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    //find the interpolation data for this t value
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

Spline::InterpolatedPT LoopingCubicHermiteSpline::getTangent(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    //find the interpolation data for this t value
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
        computePosition(t, segment),
        computeTangent(t, segment)
        );
}

Spline::InterpolatedPTC LoopingCubicHermiteSpline::getCurvature(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    //find the interpolation data for this t value
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
        computePosition(t, segment),
        computeTangent(t, segment),
        computeCurvature(t, segment)
        );
}

Spline::InterpolatedPTCW LoopingCubicHermiteSpline::getWiggle(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(t, segment)
                );
}

bool LoopingCubicHermiteSpline::isLooping(void) const
{
    return true;
}
