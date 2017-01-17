#pragma once

#include <unordered_map>
#include <vector>
#include <cmath>

namespace SplineCommon
{
    //compute the T values for the given points, with the given alpha.
    //the distance in T between adjacent points is the magitude of the distance, raised to the power alpha
    template<class InterpolationType, typename floating_t>
    floating_t computeTDiff(InterpolationType p1, InterpolationType p2, floating_t alpha);


    //compute t values for the given points, based on the alpha value
    //if innerPadding > 0, the first 'innerPadding' values will be negative, and the (innerPadding+1)'th value will be 0
    //so the spline will effectively begin at innerPadding + 1
    //this is used for splines like catmull-rom, where the first and last point are used ONLY to calculate tangent
    template<class InterpolationType, typename floating_t>
    std::unordered_map<int, floating_t> computeTValuesWithInnerPadding(
            const std::vector<InterpolationType> &points,
            floating_t alpha,
            int innerPadding
            );

    //compute compute t values for the given points, based on the alpha value
    //if outerPadding > 0, 'outerPadding' values will be added to both the beginning and ending
    //in addition to computing one t value for each point
    template<class InterpolationType, typename floating_t>
    std::unordered_map<int, floating_t> computeTValuesWithOuterPadding(
            const std::vector<InterpolationType> &points,
            floating_t alpha,
            int outerPadding
            );



    //compute the T values for the given points, with the given alpha, for use in a looping spline
    //if padding is zero, this method will return points.size() + 1 points
    //the "extra" point is because the first point in the list is represented at the beginning AND end
    //if padding is > 0, this method will also compute "extra" T values before the beginning and after the end
    //these won't actually add any extra information, but help simplify calculations that wrap around the loop
    template<class InterpolationType, typename floating_t>
    std::unordered_map<int, floating_t> computeLoopingTValues(const std::vector<InterpolationType> &points, floating_t alpha, int padding);



    //given a list of knots and a t value, return the index of the knot the t value falls within
    template<typename floating_t>
    size_t getIndexForT(const std::vector<floating_t> &knotData, floating_t t);


    template<typename floating_t>
    inline floating_t wrapGlobalT(floating_t globalT, floating_t maxT)
    {
        globalT = fmod(globalT, maxT);
        if(globalT < 0)
            return globalT + maxT;
        else
            return globalT;
    }

    //compute the arc length from a to b on the given spline
    template<template <class, typename> class Spline, class InterpolationType, typename floating_t>
    floating_t arcLength(const Spline<InterpolationType, floating_t>& spline, floating_t a, floating_t b);

    //compute the arc length from the beginning to the end on the given spline
    template<template <class, typename> class Spline, class InterpolationType, typename floating_t>
    floating_t totalLength(const Spline<InterpolationType, floating_t>& spline);
}


template<class InterpolationType, typename floating_t>
floating_t SplineCommon::computeTDiff(InterpolationType p1, InterpolationType p2, floating_t alpha)
{
    if(alpha == 0)
    {
        return 1;
    }
    else
    {
        auto distanceSq = (p1 - p2).lengthSquared();

        //if these points are right on top of each other, don't bother with the power calculation
        if(distanceSq < .0001)
        {
            return 0;
        }
        else
        {
            //multiply alpha by 0.5 so that we tke the square root of distanceSq
            //ie: result = distance ^ alpha, and distance = (distanceSq)^(0.5)
            //so: result = (distanceSq^0.5)^(alpha) = (distanceSq)^(0.5*alpha)
            //this way we don't have to do a pow AND a sqrt
            return pow(distanceSq, alpha * 0.5);
        }
    }
}

template<class InterpolationType, typename floating_t>
std::unordered_map<int, floating_t> SplineCommon::computeTValuesWithInnerPadding(
        const std::vector<InterpolationType> &points,
        floating_t alpha,
        int innerPadding
        )
{
    int size = points.size();
    int endPaddingIndex = size - 1 - innerPadding;
    int desiredMaxT = size - 2 * innerPadding - 1;

    std::unordered_map<int, floating_t> indexToT;
    indexToT.reserve(size);

    //we know points[padding] will have a t value of 0
    indexToT[innerPadding] = 0;

    //loop backwards from padding to give the earlier points negative t values
    for(int i = innerPadding - 1; i >= 0; i--)
    {
        //Points inside the padding will not be interpolated
        //so give it a negative t value, so that the first actual point can have a t value of 0
        indexToT[i] = indexToT.at(i + 1) - computeTDiff(points.at(i), points.at(i + 1), alpha);
    }

    //compute the t values of the other points
    for(int i = innerPadding + 1; i < size; i++)
    {
        indexToT[i] = indexToT.at(i - 1) + computeTDiff(points.at(i), points.at(i - 1), alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    floating_t maxTRaw = indexToT.at(endPaddingIndex);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    floating_t multiplier = desiredMaxT / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    return indexToT;
}

template<class InterpolationType, typename floating_t>
std::unordered_map<int, floating_t> SplineCommon::computeTValuesWithOuterPadding(
        const std::vector<InterpolationType> &points,
        floating_t alpha,
        int outerPadding
        )
{
    int size = points.size();

    std::unordered_map<int, floating_t> indexToT;
    indexToT.reserve(size + outerPadding * 2);

    //compute the t values each point
    indexToT[0] = 0;
    for(int i = 1; i < size; i++)
    {
        indexToT[i] = indexToT.at(i - 1) + computeTDiff(points.at(i), points.at(i - 1), alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    floating_t maxTRaw = indexToT.at(size - 1);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    floating_t desiredMaxT = size - 1;
    floating_t multiplier = desiredMaxT / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    //add padding in addition to the points - 2 on each end
    //we calculate the padding by taking the difference in t between the nearest real T values
    for(int i = size; i < outerPadding + size; i++)
    {
        floating_t tDiff = indexToT.at(i - 1) - indexToT.at(i - 2);
        indexToT[i] = indexToT.at(i - 1) + tDiff;
    }
    for(int i = -1; i >= -outerPadding; i--)
    {
        floating_t tDiff = indexToT.at(i + 2) - indexToT.at(i + 1);
        indexToT[i] = indexToT.at(i + 1) - tDiff;
    }

    return indexToT;
}

template<class InterpolationType, typename floating_t>
std::unordered_map<int, floating_t> SplineCommon::computeLoopingTValues(
        const std::vector<InterpolationType> &points,
        floating_t alpha,
        int padding)
{
    int size = points.size();

    std::unordered_map<int, floating_t> indexToT;
    indexToT.reserve(size + padding * 2);

    //compute the t values each point
    indexToT[0] = 0;
    for(int i = 1; i < size + 1; i++)
    {
        indexToT[i] = indexToT.at(i - 1) + computeTDiff(points.at(i%size), points.at(i - 1), alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    floating_t maxTRaw = indexToT.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    floating_t desiredMaxT = size;
    floating_t multiplier = desiredMaxT / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    //add padding in addition to the points - 2 on each end
    //we calculate the padding by basically wraping the difference in T values
    for(int i = 1; i < padding + 1; i++)
    {
        floating_t tDiff = indexToT.at(i) - indexToT.at(i - 1);
        indexToT[i + size] = indexToT.at(i + size - 1) + tDiff;
    }
    for(int i = -1; i >= -padding; i--)
    {
        floating_t tDiff = indexToT.at(i + size + 1) - indexToT.at(i + size);
        indexToT[i] = indexToT.at(i + 1) - tDiff;
    }

    return indexToT;
}


template<typename floating_t>
size_t SplineCommon::getIndexForT(const std::vector<floating_t> &knotData, floating_t t)
{
    //we want to find the segment whos t0 and t1 values bound x

    //if no segments bound x, return -1
    if(t <= knotData.front())
        return 0;
    if(t >= knotData.back())
        return knotData.size() - 1;

    //our initial guess will be to subtract the minimum t value, then take the floor
    int currentIndex = std::floor(t - knotData.front());
    int size = knotData.size();

    //move left or right in the array until we've found the correct index
    int searchSize = 1;
    while(t < knotData[currentIndex])
    {
        while(currentIndex >= 0 && t < knotData[currentIndex])
        {
            searchSize++;
            currentIndex -= searchSize;
        }
        if(currentIndex < 0 || t > knotData[currentIndex + 1])
        {
            currentIndex += searchSize;
            searchSize /= 4;
        }

    }
    while(t >= knotData[currentIndex + 1])
    {
        while(currentIndex < size && t >= knotData[currentIndex])
        {
            searchSize++;
            currentIndex += searchSize;
        }
        if(currentIndex >= size || t < knotData[currentIndex])
        {
            currentIndex -= searchSize;
            searchSize /= 4;
        }
    }
    return currentIndex;
}

//compute the arc length from a to b on the given spline
template<template <class, typename> class Spline, class InterpolationType, typename floating_t>
floating_t SplineCommon::arcLength(const Spline<InterpolationType, floating_t>& spline, floating_t a, floating_t b)
{
    //get the knot indices for the beginning and end
    size_t aIndex = spline.segmentForT(a);
    size_t bIndex = spline.segmentForT(b);

    //if a and b occur inside the same segment, compute the length within that segment
    //but excude cases where a > b, because that means we need to wrap around
    if(aIndex == bIndex) {
        floating_t segmentBegin = spline.segmentT(aIndex);
        floating_t segmentEnd = spline.segmentT(aIndex + 1);

        floating_t aPercent = (a - segmentBegin) / (segmentEnd - segmentBegin);
        floating_t bPercent = (b - segmentBegin) / (segmentEnd - segmentBegin);

        return spline.segmentArcLength(aIndex, aPercent, bPercent);
    }
    else {
        //a and b occur in different segments, so compute one length for every segment
        floating_t result{0};

        floating_t aBegin = spline.segmentT(aIndex);
        floating_t aEnd = spline.segmentT(aIndex + 1);
        floating_t bBegin = spline.segmentT(bIndex);
        floating_t bEnd = spline.segmentT(bIndex + 1);

        //first segment
        floating_t aPercent = (a - aBegin) / (aEnd - aBegin);
        result += spline.segmentArcLength(aIndex, aPercent, 1);

        //middle segments
        for(size_t i = aIndex + 1; i < bIndex; i++) {
            result += spline.segmentArcLength(i, 0, 1);
        }

        //last segment
        floating_t bPercent = (b - bBegin) / (bEnd - bBegin);
        result += spline.segmentArcLength(bIndex, 0, bPercent);

        return result;
    }
}

//compute the arc length from the beginning to the end on the given spline
template<template <class, typename> class Spline, class InterpolationType, typename floating_t>
floating_t SplineCommon::totalLength(const Spline<InterpolationType, floating_t>& spline)
{
    floating_t result{0};
    for(size_t i = 0; i < spline.segmentCount(); i++) {
        result += spline.segmentArcLength(i, 0, 1);
    }
    return result;
}
