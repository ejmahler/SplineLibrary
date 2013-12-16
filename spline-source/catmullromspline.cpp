
#include "catmullromspline.h"

#include <cmath>
#include <cassert>

template<class T>
T actualMod(T i, T n)
{
	int multiplier = int(i / n);
	T result = i - n * multiplier;
	return result >= 0 ? result : result + n;
}


CatmullRomSpline::CatmullRomSpline(const std::vector<Vector3D> &points, double alpha)
	:points(points)
{
    assert(points.size() >= 4);

	std::map<int, double> indexToT_Raw;
	std::map<int, Vector3D> pointMap;

	int numTotalPoints = points.size();
	int numUsedPoints;
	int firstUsedPoint;

	//if the first point is equal to the last point, this is a loop
	if(points[0] == points[points.size() - 1])
	{
		numSegments = numTotalPoints - 1;
		numUsedPoints = numTotalPoints - 1;
		firstUsedPoint = 0;
		m_isLoop = true;

		//points[-1] is a control point, so give it a negative t value, so that the first actual point can have a t value of 0
		double negativeDistance = (points.at(points.size() - 1) - points.at(points.size() - 2)).length();
		double minusOneT = -pow(negativeDistance, alpha);
		indexToT_Raw[-1] = minusOneT;
		pointMap[-1] = points[points.size() - 2];

		//we know the T value of points[0] will be 0
		indexToT_Raw[0] = 0;
		pointMap[0] = points[0];

		//compute the t values of the other points
		for(int i = 1; i < numTotalPoints + 1; i++) 
		{
			double distance = (points.at(actualMod(i, numTotalPoints)) - points.at(i - 1)).length();
			indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);

			pointMap[i] = points.at(actualMod(i, numTotalPoints));
		}
	}
	else
	{
		//this spline is not a loop
		numSegments = numTotalPoints - 3;
		numUsedPoints = numTotalPoints - 2;
		firstUsedPoint = 1;
		m_isLoop = false;

		//points[0] is a control point, so give it a negative t value, so that the first actual point can have a t value of 0
		double negativeDistance = (points.at(0) - points.at(1)).length();
		double minusOneT = -pow(negativeDistance, alpha);
		indexToT_Raw[0] = minusOneT;
		pointMap[0] = points[0];

		//compute the t values of the other points
		for(int i = 1; i < numTotalPoints; i++) 
		{
			double distance = (points.at(i) - points.at(i - 1)).length();
			indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);

			pointMap[i] = points.at(i);
		}
	}
	
	//we want to know the t value of the last segment so that we can normalize them all
	float maxTRaw = indexToT_Raw.at(firstUsedPoint + numSegments);

	//now that we have all ouf our t values and indexes figured out, normalize the t values by dividing tem by maxT
	for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
	{
		indexToT[it->first] = numSegments * it->second / maxTRaw;
	}
	maxT = indexToT.at(firstUsedPoint + numSegments);

	//compute the tangents
	std::map<int, Vector3D> tangents;
	int lastPoint = firstUsedPoint + numUsedPoints;
	for(int i = firstUsedPoint; i < lastPoint; i++)
	{
		double tPrev = indexToT.at(i - 1);
		double tCurrent = indexToT.at(i);
		double tNext = indexToT.at(i + 1);

		Vector3D pPrev = pointMap.at(i - 1);
		Vector3D pCurrent = pointMap.at(i);
		Vector3D pNext = pointMap.at(i + 1);
		
		//the tangent is the standard catmull-rom spline tangent calculation
		tangents[i] = pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
					+ pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

				 //plus a little something extra - this is derived from the pyramid contruction
				 //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0, 
				 //yielding the standard catmull-rom formula
					- pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
	}

	//pre-arrange the data needed for interpolation
	int lastSegment = firstUsedPoint + numSegments;
	for(int i = firstUsedPoint; i < lastSegment; i++) 
	{
		InterpolationData segment;
		
		segment.t0 = indexToT.at(i);
		segment.t1 = indexToT.at(i + 1);

		segment.p0 = pointMap.at(i);
		segment.p1 = pointMap.at(i + 1);

		double tDistance = segment.t1 - segment.t0;
		segment.tDistanceInverse = 1 / tDistance;

		//we scale the tangents by this segment's t distance, because wikipedia says so
		segment.m0 = tangents.at(i) * tDistance;
		segment.m1 = tangents.at(actualMod(i + 1, lastPoint)) * tDistance;

		segmentData.push_back(segment);
	}
}

CatmullRomSpline::~CatmullRomSpline()
{

}

Vector3D CatmullRomSpline::getPosition(double x) const
{
	//use modular arithmetic to bring t into an allowable range
	if(m_isLoop)
		x = actualMod(x, double(numSegments));

	//find the interpolation data for this t value
	InterpolationData segment = segmentData.at(getSegmentIndex(x));

	//we're going use the 'basis function' formula
	// http://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolation_on_an_arbitrary_interval
	double t = (x - segment.t0) * segment.tDistanceInverse;
	double oneMinusT = 1 - t;

	double basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
	double basis10 = t * oneMinusT * oneMinusT;
	double basis01 = t * t * (3 - 2*t);
	double basis11 = t * t * -oneMinusT;

	return basis00 * segment.p0 + basis10 * segment.m0 + basis01 * segment.p1 + basis11 * segment.m1;
}

InterpolatedPV CatmullRomSpline::getPositionVelocity(double x) const
{
	//use modular arithmetic to bring t into an allowable range
	if(m_isLoop)
		x = actualMod(x, double(numSegments));

	//find the interpolation data for this t value
	InterpolationData segment = segmentData.at(getSegmentIndex(x));

	//we're going use the 'basis function' formula
	// http://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolation_on_an_arbitrary_interval
	double t = (x - segment.t0) * segment.tDistanceInverse;
	double oneMinusT = 1 - t;

	double basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
	double basis10 = t * oneMinusT * oneMinusT;
	double basis01 = t * t * (3 - 2*t);
	double basis11 = t * t * -oneMinusT;

	//calculate the velocity of the spline at t. 
	//ie just compute the derivatives of all the basis functions
	double d_basis00 = 6 * t * -oneMinusT;
	double d_basis10 = (3*t - 1) * -oneMinusT;
	double d_basis01 = -d_basis00;
	double d_basis11 = t * (3 * t - 2);

	return InterpolatedPV(
		basis00 * segment.p0 +    basis10 * segment.m0 + 
		basis01 * segment.p1 +    basis11 * segment.m1,//position

		d_basis00 * segment.p0 +  d_basis10 * segment.m0 +
		d_basis01 * segment.p1 +  d_basis11 * segment.m1);//velocity
}

InterpolatedPVA CatmullRomSpline::getPositionVelocityAcceleration(double x) const
{
	//use modular arithmetic to bring t into an allowable range
	if(m_isLoop)
		x = actualMod(x, double(numSegments));

	//find the interpolation data for this t value
	InterpolationData segment = segmentData.at(getSegmentIndex(x));

	//we're going use the 'basis function' formula
	// http://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolation_on_an_arbitrary_interval
	double t = (x - segment.t0) * segment.tDistanceInverse;
	double oneMinusT = 1 - t;

	double basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
	double basis10 = t * oneMinusT * oneMinusT;
	double basis01 = t * t * (3 - 2*t);
	double basis11 = t * t * -oneMinusT;

	//calculate the velocity of the spline at t. 
	//ie just compute the derivatives of all the basis functions
	double d_basis00 = 6 * t * -oneMinusT;
	double d_basis10 = (3*t - 1) * -oneMinusT;
	double d_basis01 = -d_basis00;
	double d_basis11 = t * (3 * t - 2);

	//calculate the acceleration of the spline at t. 
	//ie just compute the second derivatives of all the basis functions
	double d2_basis00 = 6 * (2 * t - 1);
	double d2_basis10 = 2 * (3 * t - 2);
	double d2_basis01 = -d2_basis00;
	double d2_basis11 = 2 * (3 * t - 1);

	return InterpolatedPVA(
		basis00 * segment.p0 +    basis10 * segment.m0 + 
		basis01 * segment.p1 +    basis11 * segment.m1,//position

		d_basis00 * segment.p0 +  d_basis10 * segment.m0 +
		d_basis01 * segment.p1 +  d_basis11 * segment.m1,//velocity

		d2_basis00 * segment.p0 + d2_basis10 * segment.m0 +
		d2_basis01 * segment.p1 + d2_basis11 * segment.m1);//acceleration
}

double CatmullRomSpline::getT(int index) const
{
	return indexToT.at(index);
} 

double CatmullRomSpline::getMaxT(void) const
{
	return maxT;
}

int CatmullRomSpline::getNumSegments(void) const
{
	return numSegments;
}

const std::vector<Vector3D> &CatmullRomSpline::getPoints(void) const
{
	return points;
}


bool CatmullRomSpline::isLoop(void) const
{
	return m_isLoop;
}

int CatmullRomSpline::getSegmentIndex(double x) const
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
