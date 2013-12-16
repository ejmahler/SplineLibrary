#ifndef CUBICHERMITESPLINE_H
#define CUBICHERMITESPLINE_H

#include <vector>
#include <map>

#include "spline.h"

class CatmullRomSpline : public Spline
{
public:

	CatmullRomSpline(const std::vector<Vector3D> &points, double alpha = 0.5);
	~CatmullRomSpline();

	Vector3D getPosition(double x) const;
	InterpolatedPV getPositionVelocity(double x) const;
	InterpolatedPVA getPositionVelocityAcceleration(double x) const;

	double getT(int index) const;
	double getMaxT(void) const;
	int getNumSegments(void) const;

	const std::vector<Vector3D> &getPoints(void) const;

	bool isLoop(void) const;

private:
	int numSegments;
	double maxT;
	bool m_isLoop;
	
	//original point data
	std::vector<Vector3D> points;
	
	//map from index to t value. it's a map and not an array so we can store negative indexes
	std::map<int,double> indexToT;

	//a map from T1 value to all the data needed to interpolate a t value
	//there will be lots of duplication of data here, 
	//but precomputing this really speeds up the interpolation
	struct InterpolationData;
	inline int getSegmentIndex(double x) const;
	std::vector<InterpolationData> segmentData;

};

struct CatmullRomSpline::InterpolationData {
	//t values
	double t0;
	double t1;

	//points
	Vector3D p0;
	Vector3D p1;

	//tangents
	Vector3D m0;
	Vector3D m1;

	//reciprocal of distance in T between p0 and p1
	double tDistanceInverse;

	//padding to make sure we're aligned to 16 bytes
	double padding;
};

#endif // CUBICHERMITESPLINE_H
