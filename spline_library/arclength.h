#pragma once

#include <boost/math/tools/roots.hpp>

namespace ArcLength
{
    //compute the arc length from a to b on the given spline
    template<template <class, typename> class Spline, class InterpolationType, typename floating_t>
    floating_t arcLength(const Spline<InterpolationType, floating_t>& spline, floating_t a, floating_t b)
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
    floating_t totalLength(const Spline<InterpolationType, floating_t>& spline)
    {
        floating_t result{0};
        for(size_t i = 0; i < spline.segmentCount(); i++) {
            result += spline.segmentArcLength(i, 0, 1);
        }
        return result;
    }

    //compute b such that arcLength(a,b) == desiredLength
    template<template <class, typename> class Spline, class InterpolationType, typename floating_t>
    floating_t solveLength(const Spline<InterpolationType, floating_t>& spline, floating_t a, floating_t desiredLength)
    {
        size_t aIndex = spline.segmentForT(a);

        floating_t aBegin = spline.segmentT(aIndex);
        floating_t aEnd = spline.segmentT(aIndex + 1);
        floating_t aPercent = (a - aBegin) / (aEnd - aBegin);

        floating_t aLength = spline.segmentArcLength(aIndex, aPercent, 1);

        //if aLength is greater than desiredLength, then this will begin and end inside the same segment
        if(aLength > desiredLength)
        {
            auto solveFunction = [&](floating_t bPercent) {
                return spline.segmentArcLength(aIndex, aPercent, bPercent) - desiredLength;
            };

            //we can use the length from a to the segment end to formulate a pretty solid guess
            //if desired length is x% of the aLength, then our guess will be x% of the way from aPercent to 1
            floating_t desiredPercent = desiredLength / aLength;
            floating_t bGuess = aPercent + desiredPercent * (1 - aPercent);

            boost::uintmax_t max_iter = 40;
            auto result = boost::math::tools::bracket_and_solve_root(solveFunction, bGuess, floating_t(1.25), true, boost::math::tools::eps_tolerance<floating_t>(), max_iter);
            floating_t bPercent = (result.first + result.second) / 2;

            return aBegin + bPercent * (aEnd - aBegin);
        }

        //we know b is somewhere past A's segment
        desiredLength -= aLength;

        //scan through the middle segments, stopping when we reach the end of the spline or we reach the segment that contains b
        size_t bIndex = aIndex;
        floating_t bLength;
        while(++bIndex < spline.segmentCount())
        {
            bLength = spline.segmentArcLength(bIndex, 0, 1);

            if(bLength < desiredLength)
            {
                desiredLength -= bLength;
            }
            else
            {
                break;
            }
        }

        //if bIndex is equal to the segment count, we've hit the end of the spline, so return maxT
        if(bIndex == spline.segmentCount()) {
            return spline.getMaxT();
        }

        //we now know our answer lies somewhere in the segment bIndex

        //we can use the lengths we've calculated to formulate a pretty solid guess
        //if desired length is x% of the bLength, then our guess will be x% of the way from 0 to 1
        floating_t bGuess = desiredLength / bLength;

        auto solveFunction = [&](floating_t bPercent) {
            return spline.segmentArcLength(bIndex, 0, bPercent) - desiredLength;
        };

        boost::uintmax_t max_iter = 40;
        auto result = boost::math::tools::bracket_and_solve_root(solveFunction, bGuess, floating_t(1.25), true, boost::math::tools::eps_tolerance<floating_t>(), max_iter);
        floating_t bPercent = (result.first + result.second) / 2;

        floating_t bBegin = spline.segmentT(bIndex);
        floating_t bEnd = spline.segmentT(bIndex + 1);

        return bBegin + bPercent * (bEnd - bBegin);
    }
}
