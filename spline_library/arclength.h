#pragma once

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
}
