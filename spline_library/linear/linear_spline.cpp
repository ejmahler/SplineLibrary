#include "linear_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

LinearSpline::LinearSpline(const std::vector<Vector3D> &points, double alpha)
    :Spline(points)
{
    assert(points.size() >= 2);

    int size = points.size();
    numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, 0);
    maxT = indexToT.at(numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        LinearSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D LinearSpline::getPosition(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

Spline::InterpolatedPT LinearSpline::getTangent(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent()
                );
}

Spline::InterpolatedPTC LinearSpline::getCurvature(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTC(
                segment.computePosition(localT),
                segment.computeTangent(),
                Vector3D() //curvature is always 0 for linear spline
                );
}

Spline::InterpolatedPTCW LinearSpline::getWiggle(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTCW(
                segment.computePosition(localT),
                segment.computeTangent(),
                Vector3D(), //curvature and wiggle are always 0 for linear spline
                Vector3D()
                );
}

double LinearSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LinearSpline::getMaxT(void) const
{
    return maxT;
}

bool LinearSpline::isLooping(void) const
{
    return false;
}
