#include "cubic_hermite_spline.h"

Vector3D CubicHermiteSpline::getPosition(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

InterpolatedPT CubicHermiteSpline::getTangent(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

InterpolatedPTC CubicHermiteSpline::getCurvature(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

int CubicHermiteSpline::getSegmentIndex(double x) const
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

double CubicHermiteSpline::getT(int index) const
{
    return indexToT.at(index);
}

double CubicHermiteSpline::getMaxT(void) const
{
    return maxT;
}

int CubicHermiteSpline::getNumSegments(void) const
{
    return numSegments;
}

const std::vector<Vector3D> &CubicHermiteSpline::getPoints(void) const
{
    return points;
}

bool CubicHermiteSpline::isLooping(void) const
{
    return false;
}
