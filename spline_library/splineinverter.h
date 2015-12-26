#ifndef SplineInverter_H
#define SplineInverter_H

#include <vector>
#include <array>
#include <memory>

#include <boost/math/tools/minima.hpp>

#include "spline_library/spline.h"
#include "spline_library/utils/splinesample_adaptor.h"

template<class InterpolationType, typename floating_t, int sampleDimension>
std::array<floating_t, sampleDimension> convertPoint(const InterpolationType& point);

template<class InterpolationType, typename floating_t=float, int sampleDimension=2>
class SplineInverter
{
public:
    SplineInverter(const Spline<InterpolationType, floating_t> &spline, int samplesPerT = 10);

    floating_t findClosestT(const InterpolationType &queryPoint) const;

private: //methods
    SplineSamples<sampleDimension, floating_t> makeSplineSamples(void) const;

private: //data
    const Spline<InterpolationType, floating_t> &spline;

    //distance in t between samples
    floating_t sampleStep;

    SplineSampleTree<sampleDimension, floating_t> sampleTree;
};

template<class InterpolationType, typename floating_t, int sampleDimension>
SplineInverter<InterpolationType, floating_t, sampleDimension>::SplineInverter(
        const Spline<InterpolationType, floating_t> &spline,
        int samplesPerT)
    :spline(spline),
      sampleStep(1.0 / floating_t(samplesPerT)),
      sampleTree(makeSplineSamples())
{

}

template<class InterpolationType, typename floating_t, int sampleDimension>
SplineSamples<sampleDimension, floating_t> SplineInverter<InterpolationType, floating_t, sampleDimension>::makeSplineSamples(void) const
{
    SplineSamples<sampleDimension, floating_t> samples;

    //first step is to populate the splineSamples map
    //we're going to have sampled T values sorted by x coordinate
    floating_t currentT = 0;
    floating_t maxT = spline.getMaxT();
    while(currentT < maxT)
    {
        auto sampledPoint = convertPoint<InterpolationType, floating_t, sampleDimension>(spline.getPosition(currentT));

        samples.pts.emplace_back(sampledPoint, currentT);

        currentT += sampleStep;
    }

    //if the spline isn't a loop and the final t value isn't very very close to maxT, we have to add a sample for maxT
    floating_t lastT = samples.pts.at(samples.pts.size() - 1).t;
    if(!spline.isLooping() && std::abs(lastT / maxT - 1) > .0001)
    {
        auto sampledPoint = convertPoint<InterpolationType, floating_t, sampleDimension>(spline.getPosition(currentT));

        samples.pts.emplace_back(sampledPoint, currentT);
    }

    return samples;
}

template<class InterpolationType, typename floating_t, int sampleDimension>
floating_t SplineInverter<InterpolationType, floating_t, sampleDimension>::findClosestT(const InterpolationType &queryPoint) const
{
    auto convertedQueryPoint = convertPoint<InterpolationType, floating_t, sampleDimension>(queryPoint);
    floating_t closestSampleT = sampleTree.findClosestSample(convertedQueryPoint);

    //compute the first derivative of distance to spline at the sample point
    auto sampleResult = spline.getTangent(closestSampleT);
    InterpolationType sampleDisplacement = sampleResult.position - queryPoint;
    floating_t sampleDistanceSlope = InterpolationType::dotProduct(sampleDisplacement.normalized(), sampleResult.tangent);

    //if the spline is not a loop there are a few special cases to account for
    if(!spline.isLooping())
    {
        //if closest sample T is 0, we are on an end. so if the slope is positive, we have to just return the end
        if(closestSampleT == 0 && sampleDistanceSlope > 0)
            return 0;

        //if the closest sample T is max T we are on an end. so if the slope is negative, just return the end
        if(closestSampleT == spline.getMaxT() && sampleDistanceSlope < 0)
            return spline.getMaxT();
    }

    //step forwards or backwards in the spline until we find a point where the distance slope has flipped sign.
    //because "currentsample" is the closest point, the "next" sample's slope MUST have a different sign
    //otherwise that sample would be closer
    //note: this assumption is only true if the samples are close together

    //if sample distance slope is positive we want to move backwards in t, otherwise forwards
    floating_t a, b;
    if(sampleDistanceSlope > 0)
    {
        a = closestSampleT - sampleStep;
        b = closestSampleT;
    }
    else
    {
        a = closestSampleT;
        b = closestSampleT + sampleStep;
    }

    auto distanceFunction = [&spline = spline, queryPoint](floating_t t) {
        return (spline.getPosition(t) - queryPoint).lengthSquared();
    };

    //we know that the actual closest T is now between a and b
    //use brent's method to find the actual closest point, using a and b as bounds
    auto result = boost::math::tools::brent_find_minima(distanceFunction, a, b, 16);
    return result.first;
}

#endif // SplineInverter_H
