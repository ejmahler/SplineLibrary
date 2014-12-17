#ifndef QUINTICHERMITESPLINE_H
#define QUINTICHERMITESPLINE_H

#include <unordered_map>

#include "spline_library/spline.h"
#include "spline_library/hermite/quintic/quintic_hermite_spline_kernel.h"

class QuinticHermiteSpline final : public Spline
{
//constructors
public:
    QuinticHermiteSpline(const std::vector<Vector3D> &points,
                         const std::vector<Vector3D> &tangents,
                         const std::vector<Vector3D> &curvatures,
                         float alpha = 0.0f
                         );
    QuinticHermiteSpline(const std::vector<Vector3D> &points, double alpha = 0.0f);

//methods
public:
    Vector3D getPosition(double x) const override;
    InterpolatedPT getTangent(double x) const override;
    InterpolatedPTC getCurvature(double x) const override;
    InterpolatedPTCW getWiggle(double x) const override;

    double getT(int index) const override;
    double getMaxT(void) const override;

    bool isLooping(void) const override;

//data
protected:
	//a vector containing pre-computed datasets, one per segment
    //there will be lots of duplication of data here,
    //but precomputing this really speeds up the interpolation
	int numSegments;
    std::vector<QuinticHermiteSplineKernel::InterpolationData<Vector3D>> segmentData;

    double maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,double> indexToT;
};

#endif // QUINTICHERMITESPLINE_H
