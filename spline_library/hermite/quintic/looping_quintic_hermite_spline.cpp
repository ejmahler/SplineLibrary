#include "looping_quintic_hermite_spline.h"

#include "spline_library/utils/t_calculator.h"

#include <cmath>
#include <cassert>

LoopingQuinticHermiteSpline::LoopingQuinticHermiteSpline()
{
}

LoopingQuinticHermiteSpline::LoopingQuinticHermiteSpline(
        const std::vector<Vector3D> &points,
        const std::vector<Vector3D> &tangents,
        const std::vector<Vector3D> &curvatures,
        float alpha
        )
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());
    assert(points.size() == curvatures.size());

    this->points = points;

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 0;
    indexToT = TCalculator::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangents.at(i) * tDistance;
        segment.m1 = tangents.at((i + 1)%size) * tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.c0 = curvatures.at(i) * tDistance * tDistance * tDistance;
        segment.c1 = curvatures.at((i + 1)%size) * tDistance * tDistance * tDistance;

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
