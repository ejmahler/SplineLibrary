#include "looping_quintic_hermite_spline.h"

#include <cmath>
#include <cassert>

LoopingQuinticHermiteSpline::LoopingQuinticHermiteSpline()
{
}

LoopingQuinticHermiteSpline::LoopingQuinticHermiteSpline(
        const std::vector<Vector3D> &points,
        const std::vector<Vector3D> &tangents,
        const std::vector<Vector3D> &curvatures
        )
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());
    assert(points.size() == curvatures.size());

    this->points = points;

    //i would love to be able to support changing alphas for quintic catmull rom splines!
    //but there's no literature whatsoever on how to choose tangents when t values are unevenly spaced
    //if you know how to do it, let me know or make a pull request :D
    //until then we're just going to hardcode alpha to 0 and make sure nothing in the code
    //besides tangent selection assumes t values are evenly spaced
    float alpha = 0;

    std::unordered_map<int, double> indexToT_Raw;

    int size = points.size();
    numSegments = size;

    //we know points[0] will have a t value of 0
    indexToT_Raw[0] = 0;

    //compute the t values of the other points
    for(int i = 1; i < size + 1; i++)
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

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.c0 = curvatures.at(i) * tDistance * tDistance * tDistance;
        segment.c1 = curvatures.at(i + 1) * tDistance * tDistance * tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D LoopingQuinticHermiteSpline::getPosition(double x) const
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

Spline::InterpolatedPT LoopingQuinticHermiteSpline::getTangent(double x) const
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

Spline::InterpolatedPTC LoopingQuinticHermiteSpline::getCurvature(double x) const
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

Spline::InterpolatedPTCW LoopingQuinticHermiteSpline::getWiggle(double x) const
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

bool LoopingQuinticHermiteSpline::isLooping(void) const
{
    return true;
}
