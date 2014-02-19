#include "linear_spline.h"

#include "spline_library/utils/t_calculator.h"

#include <cmath>
#include <cassert>

LinearSpline::LinearSpline()
{

}

LinearSpline::LinearSpline(const std::vector<Vector3D> &points, double alpha)
{
    assert(points.size() >= 2);

    this->points = points;

    int size = points.size();
    numSegments = size - 1;

    //compute the T values for each point
    indexToT = TCalculator::computeTValues(points, alpha, 0);
    maxT = indexToT.at(numSegments);

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

        segmentData.push_back(segment);
    }
}

Vector3D LinearSpline::getPosition(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

Spline::InterpolatedPT LinearSpline::getTangent(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(segment)
                );
}

Spline::InterpolatedPTC LinearSpline::getCurvature(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(segment),
                Vector3D() //curvature is always 0 for linear spline
                );
}

Spline::InterpolatedPTCW LinearSpline::getWiggle(double x) const
{
    int index = getSegmentIndex(x);
    InterpolationData segment = segmentData.at(index);
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(segment),
                Vector3D(), //curvature and wiggle are always 0 for linear spline
                Vector3D()
                );
}

int LinearSpline::getSegmentIndex(double x) const
{
    //we want to find the segment whos t0 and t1 values bound x

    //if no segments bound x, return -1
    if(x < segmentData[0].t0)
        return 0;
    if(x > segmentData[numSegments - 1].t1)
        return numSegments - 1;

    //perform a binary search on segmentData
    int currentMin = 0;
    int currentMax = segmentData.size() - 1;
    int currentIndex = (currentMin + currentMax) / 2;

    //keep looping as long as this segment does not bound x

    while(x < segmentData[currentIndex].t0 || x > segmentData[currentIndex].t1)
    {
        //if t0 is greater than x, search the left half of the array
        if(segmentData[currentIndex].t0 > x)
        {
            currentMax = currentIndex - 1;
        }

        //the only other possibility is that t1 is less than x, so search the right half of the array
        else
        {
            currentMin = currentIndex + 1;
        }
        currentIndex = (currentMin + currentMax) / 2;
    }
    return currentIndex;
}

double LinearSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LinearSpline::getMaxT(void) const
{
    return maxT;
}

const std::vector<Vector3D> &LinearSpline::getPoints(void) const
{
    return points;
}

bool LinearSpline::isLooping(void) const
{
    return false;
}
