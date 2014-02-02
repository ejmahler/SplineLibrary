#include "cubic_b_spline.h"

#include "../utils/t_calculator.h"

#include <cmath>
#include <cassert>

CubicBSpline::CubicBSpline()
{

}

CubicBSpline::CubicBSpline(const std::vector<Vector3D> &points)
{
    assert(points.size() >= 4);
    double alpha = 0.0;

    this->points = points;

    int size = points.size();
    int padding = 1;
    numSegments = size - 3;

    //compute the T values for each point
    indexToT = TCalculator::computeTValues(points, alpha, padding);
    maxT = indexToT.at(padding + numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = padding; i < padding + numSegments; i++)
    {
        InterpolationData segment;

        segment.t1 = indexToT.at(i);
        segment.t2 = indexToT.at(i + 1);

        segment.p0 = points.at(i - 1);
        segment.p1 = points.at(i);
        segment.p2 = points.at(i + 1);
        segment.p3 = points.at(i + 2);

        segmentData.push_back(segment);
    }
}

Vector3D CubicBSpline::getPosition(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return computePosition(t, segment);
}

Spline::InterpolatedPT CubicBSpline::getTangent(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC CubicBSpline::getCurvature(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW CubicBSpline::getWiggle(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(t, segment)
                );
}

int CubicBSpline::getSegmentIndex(double x) const
{
    //we want to find the segment whos t0 and t1 values bound x

    //if no segments bound x, return -1
    if(x < segmentData[0].t1)
        return 0;
    if(x > segmentData[numSegments - 1].t2)
        return numSegments - 1;

    //perform a binary search on segmentData
    int currentMin = 0;
    int currentMax = segmentData.size() - 1;
    int currentIndex = (currentMin + currentMax) / 2;

    //keep looping as long as this segment does not bound x

    while(x < segmentData[currentIndex].t1 || x > segmentData[currentIndex].t2)
    {
        //if t0 is greater than x, search the left half of the array
        if(segmentData[currentIndex].t1 > x)
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

double CubicBSpline::getT(int index) const
{
    return indexToT.at(index);
}

double CubicBSpline::getMaxT(void) const
{
    return maxT;
}

int CubicBSpline::getNumSegments(void) const
{
    return numSegments;
}

const std::vector<Vector3D> &CubicBSpline::getPoints(void) const
{
    return points;
}

bool CubicBSpline::isLooping(void) const
{
    return false;
}
