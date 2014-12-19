#ifndef SplineInverter_H
#define SplineInverter_H

#include <vector>
#include <memory>

#include "spline_library/spline.h"
#include "spline_library/utils/splinesample_adaptor.h"
#include "utils/optimization.h"

template<class InterpolationType, typename floating_t=float>
class SplineInverter
{
public:
    SplineInverter(const std::shared_ptr<Spline<InterpolationType, floating_t>> &spline, int samplesPerT = 10);
	~SplineInverter();

    floating_t findClosestT(const InterpolationType &queryPoint) const;

private: //data
    std::shared_ptr<Spline<InterpolationType, floating_t>> spline;

    //distance in t between samples
    floating_t sampleStep;

	//when we find a T whose abs(distance slope) is less than this tolerance, we return
    floating_t slopeTolerance;

    //inner class used to provide an abstraction between the spline inverter and nanoflann
    std::unique_ptr<SplineSampleTree> sampleTree;
};

template<class InterpolationType, typename floating_t>
SplineInverter<InterpolationType, floating_t>::SplineInverter(
        const std::shared_ptr<Spline<InterpolationType,floating_t>> &spline,
        int samplesPerT)
    :spline(spline), sampleStep(1.0 / floating_t(samplesPerT)), slopeTolerance(0.01)
{
    SplineSamples3D samples;

    //first step is to populate the splineSamples map
    //we're going to have sampled T values sorted by x coordinate
    floating_t currentT = 0;
    floating_t maxT = spline->getMaxT();
    while(currentT < maxT)
    {
        auto sampledPoint = spline->getPosition(currentT);

        samples.pts.emplace_back(sampledPoint.x(), sampledPoint.y(), sampledPoint.z(), currentT);

        currentT += sampleStep;
    }

    //if the spline isn't a loop and the final t value isn't very very close to maxT, we have to add a sample for maxT
    floating_t lastT = samples.pts.at(samples.pts.size() - 1).t;
    if(!spline->isLooping() && abs(lastT / maxT - 1) > .0001)
    {
        auto sampledPoint = spline->getPosition(maxT);

        samples.pts.emplace_back(sampledPoint.x(), sampledPoint.y(), sampledPoint.z(), currentT);
    }

    //populate the sample kd-tree
    sampleTree = std::unique_ptr<SplineSampleTree>(new SplineSampleTree(samples));
}

template<class InterpolationType, typename floating_t>
SplineInverter<InterpolationType, floating_t>::~SplineInverter()
{

}

template<class InterpolationType, typename floating_t>
floating_t SplineInverter<InterpolationType, floating_t>::findClosestT(const InterpolationType &queryPoint) const
{
    floating_t closestSampleT = sampleTree->findClosestSample(queryPoint);

    //define a lambda to compute the slope of the distance to the querypoint at T
    auto splineInstance = spline;
    auto distanceSlopeFunction = [splineInstance, queryPoint](floating_t t) {
        auto result = splineInstance->getTangent(t);

        //get the displacement from the spline at T to the query point
        InterpolationType displacement = result.position - queryPoint;

        //find projection of spline velocity onto displacement
        return InterpolationType::dotProduct(displacement.normalized(), result.tangent);
    };

    floating_t sampleDistanceSlope = distanceSlopeFunction(closestSampleT);

    //if the slope is very close to 0, just return the sampled point
    if(abs(sampleDistanceSlope) < slopeTolerance)
        return closestSampleT;

    //if the spline is not a loop there are a few special cases to account for
    if(!spline->isLooping())
    {
        //if closest sample T is 0, we are on an end. so if the slope is positive, we have to just return the end
        if(abs(closestSampleT) < .0001 && sampleDistanceSlope > 0)
            return closestSampleT;

        //if the closest sample T is max T we are on an end. so if the slope is negative, just return the end
        if(abs(closestSampleT / spline->getMaxT() - 1) < .0001 && sampleDistanceSlope < 0)
            return closestSampleT;
    }

    //step forwards or backwards in the spline until we find a point where the distance slope has flipped sign.
    //because "currentsample" is the closest point, the "next" sample's slope MUST have a different sign
    //otherwise that sample would be closer
    //note: this assumption is only true if the samples are close together

    floating_t a = closestSampleT;

    //if sample distance slope is positive we want to move backwards in t, otherwise forwards
    floating_t b = closestSampleT - sampleStep * sign(sampleDistanceSlope);

    floating_t aValue = sampleDistanceSlope;
    floating_t bValue = distanceSlopeFunction(b);

    //we know that the actual closest T is now between a and b
    //use brent's method to find the actual closest point, using a and b as bounds
    return Optimization::brentsMethod(distanceSlopeFunction, a, aValue, b, bValue);
}

#endif // SplineInverter_H
