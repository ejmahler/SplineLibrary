#ifndef SPLINELENGTHCALCULATOR_H
#define SPLINELENGTHCALCULATOR_H

#include <memory>
#include <atomic>

#include "spline_library/spline.h"
class InterpolationType;

template<class InterpolationType, typename floating_t=float>
class SplineLengthCalculator
{
public:
    SplineLengthCalculator(const Spline<InterpolationType,floating_t> &spline);

    //precise but possibly slower method to find an approximation of the spline length from the begin T to the end T
    //if the spline is a looping spline and useShortestPath is true, this will try to find a shorter path around the spline by going "backwards" if possible
    //useShortestPath has no effect on non-looping splines
    floating_t findLength(floating_t beginT, floating_t endT, bool useShortestPath=false, floating_t eps = 0.0025) const;

private: //methods
    floating_t computeLength(floating_t beginT, floating_t endT, floating_t eps) const;

    //recursively compute the length of the given segment
    floating_t computeLengthHelper(floating_t beginT, const InterpolationType &beginPosition,
                                floating_t endT, const InterpolationType &endPosition, floating_t eps) const;

private: //data

    const Spline<InterpolationType,floating_t> &spline;
    floating_t maxT;
    mutable std::atomic<floating_t> atomic_splineLength;

    const static floating_t RECURSIVE_MINIMUM_INTERVAL;
    const static floating_t LENGTH_INTERVAL;
};

template<class InterpolationType, typename floating_t>
const floating_t SplineLengthCalculator<InterpolationType,floating_t>::RECURSIVE_MINIMUM_INTERVAL = 0.0001;

template<class InterpolationType, typename floating_t>
const floating_t SplineLengthCalculator<InterpolationType,floating_t>::LENGTH_INTERVAL = 0.25;

template<class InterpolationType, typename floating_t>
SplineLengthCalculator<InterpolationType,floating_t>::SplineLengthCalculator(
        const Spline<InterpolationType,floating_t> &spline)
    :
      spline(spline),
      maxT(spline.getMaxT()),

      //the spline length will be lazy-computed when it's needed
      atomic_splineLength(-1)
{

}

template<class InterpolationType, typename floating_t>
floating_t SplineLengthCalculator<InterpolationType,floating_t>::findLength(floating_t beginT, floating_t endT, bool useShortestPath, floating_t eps) const
{

    //if this is a looping spline and the calles had requested the shortest path, the behavior will be slightly different
    if(useShortestPath && spline.isLooping())
    {
        //use fmod on both args to make sure we're not going around the loop multiple times
        beginT = fmod(beginT, maxT);
        endT = fmod(endT, maxT);

        //make sure they're both positive
        beginT = (beginT < 0) ? (beginT + maxT) : beginT;
        endT = (endT < 0) ? (endT + maxT) : endT;

        //make sure beginT is less than endT
        if(beginT > endT) {
            std::swap(beginT, endT);
        }

        //compute length
        floating_t computedLength = computeLength(beginT, endT, eps);

        floating_t splineLength = atomic_splineLength.load(std::memory_order_consume);

        //lazy compute the spline length
        if(splineLength < 0)
        {
            //modifying a variable inside of a const function can be non thread safe, but in this case thread safety is not a concern
            //we're using std::atomic so the worst possible thing that could happen is this variable could be computed multiple times by simultaneous threads
            //but if that happens, it'll be the same both times so who cares
            splineLength = computeLength(0, maxT, eps);
            atomic_splineLength.store(splineLength, std::memory_order_release);
        }

        //that was one direction around, the other direction will be (splineLength - computedLength). return the smaller of the two
        return std::min(computedLength, splineLength - computedLength);
    }
    else
    {
        //make sure beginT is less than endT
        if(beginT > endT) {
            std::swap(beginT, endT);
        }

        return computeLength(beginT, endT, eps);
    }
}

template<class InterpolationType, typename floating_t>
floating_t SplineLengthCalculator<InterpolationType,floating_t>::computeLength(floating_t beginT, floating_t endT, floating_t eps) const
{
    floating_t tDistance = endT - beginT;

    //find the number of segments we're going to use
    int numSegments = ceil(tDistance / LENGTH_INTERVAL);

    //now that we know the number of segments, find the actual interval
    floating_t interval = tDistance / numSegments;

    //for each segment, compute the length of a circle arc passing though that segment's endpoints
    auto previousData = spline.getPosition(beginT);
    floating_t lengthSum = 0;

    for(int i = 0; i < numSegments; i++)
    {
        floating_t currentT = beginT + (i + 1) * interval;

        auto currentData = spline.getPosition(currentT);
        lengthSum += computeLengthHelper(
                    currentT - interval,    previousData,
                    currentT,               currentData, eps);

        previousData = currentData;
    }

    return lengthSum;
}

template<class InterpolationType, typename floating_t>
floating_t SplineLengthCalculator<InterpolationType,floating_t>::computeLengthHelper(
        floating_t beginT, const InterpolationType &beginPosition,
        floating_t endT, const InterpolationType &endPosition, floating_t eps) const
{
    //compute the midpoint between the start and end
    floating_t midT = (beginT + endT) * 0.5;
    InterpolationType midPosition = spline.getPosition(midT);

    //compute the length squared for start->mid plus mid->end
    floating_t fullLength = (endPosition - beginPosition).lengthSquared();
    floating_t halfLengths = ((endPosition - midPosition).lengthSquared() + (midPosition - beginPosition).lengthSquared()) * 2;

    //if the difference bween the total lengh squared and the length squared for the two halves is small enough, return the length of the whole segment
    //eps is a percentage;
    floating_t limit = 1 - eps;
    if(fullLength / halfLengths > limit || (endT - beginT < RECURSIVE_MINIMUM_INTERVAL))
    {
        return sqrt(fullLength);
    }
    else
    {
        //there is too much of a difference, so recursively compute the lengths
        floating_t firstHalf = computeLengthHelper(beginT, beginPosition, midT, midPosition, eps);
        floating_t secondHalf = computeLengthHelper(midT, midPosition, endT, endPosition, eps);

        return firstHalf + secondHalf;
    }
}

#endif // SPLINELENGTHCALCULATOR_H
