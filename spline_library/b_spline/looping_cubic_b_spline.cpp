#include "looping_cubic_b_spline.h"

#include <cassert>
#include <cmath>

LoopingCubicBSpline::LoopingCubicBSpline(const std::vector<Vector3D> &points)
{
    assert(points.size() >= 3);
    double alpha = 0.0;

    this->points = points;

    std::unordered_map<int, double> indexToT_Raw;
    std::unordered_map<int, Vector3D> pointMap;

    int size = points.size();
    numSegments = size;

    //we know points[0] will have a t value of 0
    indexToT_Raw[0] = 0;
    pointMap[0] = points[0];

    //loop backwards from 0 to give the earlier points negative t values
    for(int i = -1; i >= -1; i--)
    {
        //points[1] is a control point, so give it a nagative t value, so that the first actual point can have a t value of 0
        double distance = (points.at((i + size)%size) - points.at((i + 1 + size)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i + 1] - pow(distance, alpha);
        pointMap[i] = points.at((i + size)%size);
    }

    //compute the t values of the other points
    for(int i = 1; i < size + 3; i++)
    {
        double distance = (points.at(i%size) - points.at((i - 1)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);

        pointMap[i] = points.at((i + size)%size);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(size);

    //now that we have all of our t values and indexes figured out, normalize the t values by dividing tem by maxT
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = numSegments * it->second / maxTRaw;
    }
    maxT = indexToT.at(size);


    //pre-arrange the data needed for interpolation
    for(int i = 0; i < size + 1; i++)
    {
        InterpolationData segment;

        segment.t1 = indexToT.at(i);
        segment.t2 = indexToT.at(i + 1);

        segment.p0 = pointMap.at(i - 1);
        segment.p1 = pointMap.at(i);
        segment.p2 = pointMap.at(i + 1);
        segment.p3 = pointMap.at(i + 2);

        segmentData.push_back(segment);
    }
}

Vector3D LoopingCubicBSpline::getPosition(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return computePosition(t, segment);
}

Spline::InterpolatedPT LoopingCubicBSpline::getTangent(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC LoopingCubicBSpline::getCurvature(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW LoopingCubicBSpline::getWiggle(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(t, segment)
                );
}

bool LoopingCubicBSpline::isLooping(void) const
{
    return true;
}
