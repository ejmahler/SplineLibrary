#include "looping_quintic_cr_spline.h"


#include <cmath>
#include <cassert>

LoopingQuinticCRSpline::LoopingQuinticCRSpline(const std::vector<Vector3D> &points)
    :points(points)
{
    assert(points.size() >= 6);

	//i would love to be able to support changing alphas for quintic catmull rom splines!
	//but there's no literature whatsoever on how to choose tangents when t values are unevenly spaced
	//if you know how to do it, let me know or make a pull request :D
	//until then we're just going to hardcode alpha to 0 and make sure nothing in the code
    //besides tangent selection assumes t values are evenly spaced
    float alpha = 0;

    std::unordered_map<int, double> indexToT_Raw;
    std::unordered_map<int, Vector3D> pointMap;

    int size = points.size();
    numSegments = size;

	//we know points[0] will have a t value of 0
	indexToT_Raw[0] = 0;
    pointMap[0] = points[0];

	//loop backwards from 0 to give the earlier points negative t values
    for(int i = -1; i >= -3; i--)
	{
		//points[1] is a control point, so give it a nagative t value, so that the first actual point can have a t value of 0
		double distance = (points.at((i + size)%size) - points.at((i + 1 + size)%size)).length();
		indexToT_Raw[i] = indexToT_Raw[i + 1] - pow(distance, alpha);
		pointMap[i] = points.at((i + size)%size);
	}

    //compute the t values of the other points
    for(int i = 1; i < size + 4; i++)
    {
        double distance = (points.at(i%size) - points.at((i - 1)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);

        pointMap[i] = points.at((i + size)%size);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing tem by maxT
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = numSegments * it->second / maxTRaw;
    }
    maxT = indexToT.at(size);

    //compute the tangents
    std::map<int, Vector3D> tangentMap;
    for(int i = -1; i < size + 2; i++)
    {
        Vector3D pPrev = pointMap.at(i - 1);
        Vector3D pNext = pointMap.at(i + 1);

        //the tangent is the standard catmull-rom spline tangent calculation
        tangentMap[i] = (pNext - pPrev) * 0.5;
    }

    //compute the curvatures
    std::map<int, Vector3D> curveMap;
    for(int i = 0; i < size + 1; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tNext = indexToT.at(i + 1);

        Vector3D tangentPrev = tangentMap.at(i - 1);
        Vector3D tangentNext = tangentMap.at(i + 1);

        //the curve uses the same formula as the tangent
        curveMap[i] = (tangentNext - tangentPrev) / (tNext - tPrev);
	}


    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = pointMap.at(i);
        segment.p1 = pointMap.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangentMap.at(i) * tDistance;
        segment.m1 = tangentMap.at(i + 1) * tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.c0 = curveMap.at(i) * tDistance;
        segment.c1 = curveMap.at(i + 1) * tDistance * tDistance * tDistance;

        segmentData.push_back(segment);
    }
}

LoopingQuinticCRSpline::~LoopingQuinticCRSpline()
{

}

Vector3D LoopingQuinticCRSpline::getPosition(double x) const
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

InterpolatedPT LoopingQuinticCRSpline::getTangent(double x) const
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

InterpolatedPTC LoopingQuinticCRSpline::getCurvature(double x) const
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

double LoopingQuinticCRSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LoopingQuinticCRSpline::getMaxT(void) const
{
    return maxT;
}

int LoopingQuinticCRSpline::getNumSegments(void) const
{
    return numSegments;
}

const std::vector<Vector3D> &LoopingQuinticCRSpline::getPoints(void) const
{
    return points;
}


bool LoopingQuinticCRSpline::isLoop(void) const
{
    return true;
}
