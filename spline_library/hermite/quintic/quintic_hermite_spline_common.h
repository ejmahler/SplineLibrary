#pragma once

#include <vector>

#include "../../utils/spline_setup.h"
#include "../../utils/calculus.h"

template<class InterpolationType, typename floating_t>
class QuinticHermiteSplineCommon
{
public:
    struct alignas(16) QuinticHermiteSplinePoint
    {
        InterpolationType position, tangent, curvature;
    };

    inline QuinticHermiteSplineCommon(void) = default;
    inline QuinticHermiteSplineCommon(std::vector<QuinticHermiteSplinePoint> points, std::vector<floating_t> knots)
        :points(std::move(points)), knots(std::move(knots))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return computePosition(knotIndex, tDiff, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(knotIndex, tDiff, localT),
                    computeTangent(knotIndex, tDiff, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computePosition(knotIndex, tDiff, localT),
                    computeTangent(knotIndex, tDiff, localT),
                    computeCurvature(knotIndex, tDiff, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computePosition(knotIndex, tDiff, localT),
                    computeTangent(knotIndex, tDiff, localT),
                    computeCurvature(knotIndex, tDiff, localT),
                    computeWiggle(knotIndex, tDiff, localT)
                    );
    }

    inline floating_t getLength(floating_t a, floating_t b) const
    {
        //get the knot indices for the beginning and end
        size_t aIndex = size_t(a);
        size_t bIndex = size_t(b);

        size_t numSegments = knots.size() - 1;

        if(aIndex > numSegments)
            aIndex = numSegments;
        if(bIndex > numSegments)
            bIndex = numSegments;

        //if a and b occur inside the same segment, compute the length within that segment
        //but excude cases where a > b, because that means we need to wrap around
        if(aIndex == bIndex && a <= b) {
            return computeSegmentLength(aIndex, a - aIndex, b - aIndex);
        }
        else {
            //a and b occur in different segments, so compute one length for every segment
            floating_t result{0};

            //first segment
            result += computeSegmentLength(aIndex, a - aIndex, 1);

            //last segment
            result += computeSegmentLength(bIndex, 0, b - bIndex);

            //if b index is less than a index, that means the user wants to wrap around the end of the spline and back to the beginning
            //if so, add the number of points in the spline to bIndex, and we'll use mod to make sure it stays in range
            if(bIndex <= aIndex)
                bIndex += numSegments;

            //middle segments
            for(size_t i = aIndex + 1; i < bIndex; i++) {
                result += computeSegmentLength(i%numSegments, 0, 1);
            }

            return result;
        }
    }

    inline floating_t getTotalLength(void) const
    {
        floating_t result{0};
        for(size_t i = 0; i < knots.size() - 1; i++) {
            result += computeSegmentLength(i, 0, 1);
        }
        return result;
    }



private: //methods
    inline InterpolationType computePosition(size_t index, floating_t tDiff, floating_t t) const
    {
        auto oneMinusT = 1 - t;

        //this is a logical extension of the cubic hermite spline's basis functions
        //that has one basis function for t0 position, one for t1 position
        //one for t0 tangent (1st derivative of position), and one for t1 tangent
        //this adds 2 more basis functions, one for t0 curvature (2nd derivative) and t1 curvature
        //see this paper for details http://www.rose-hulman.edu/~finn/CCLI/Notes/day09.pdf
        auto basis00 = (oneMinusT * oneMinusT * oneMinusT * (t * (6 * t + 3) + 1));
        auto basis10 = (t * oneMinusT * oneMinusT * oneMinusT * (3 * t + 1));
        auto basis20 = 0.5 * oneMinusT * oneMinusT * oneMinusT * t * t;
        auto basis21 = 0.5 * oneMinusT * oneMinusT * t * t * t;
        auto basis11 = t * t * t * oneMinusT * (t * 3 - 4);
        auto basis01 = t * t * t * (t * (6 * t - 15) + 10);

        return
                basis00 * points[index].position +
                basis10 * tDiff * points[index].tangent +
                basis20 * tDiff * tDiff * points[index].curvature +

                basis21 * tDiff * tDiff * points[index + 1].curvature +
                basis11 * tDiff * points[index + 1].tangent +
                basis01 * points[index + 1].position;
    }

    inline InterpolationType computeTangent(size_t index, floating_t tDiff, floating_t t) const
    {
        auto oneMinusT = 1 - t;

        //we're computing the derivative of the computePosition function with respect to t
        //we can do this by computing the derivatives of each of its basis functions.
        //thankfully this can easily be done analytically since they're polynomials!
        auto d_basis00 = 30 * oneMinusT * (t - 1) * t * t;
        auto d_basis10 = oneMinusT * oneMinusT * (1 - 3 * t) * (5 * t + 1);
        auto d_basis20 = 0.5 * oneMinusT * oneMinusT * t * (2 - 5 * t);
        auto d_basis21 = 0.5 * oneMinusT * t * t * (3 - 5 * t);
        auto d_basis11 = -t * t * (2 - 3 * t) * (5 * t - 6);
        auto d_basis01 = 30 * oneMinusT * oneMinusT * t * t;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the derivative of the position function and nothing else
        //if you know why please let me know
        return (
                    d_basis00 * points[index].position +
                    d_basis10 * tDiff * points[index].tangent +
                    d_basis20 * tDiff * tDiff * points[index].curvature +

                    d_basis21 * tDiff * tDiff * points[index + 1].curvature +
                    d_basis11 * tDiff * points[index + 1].tangent +
                    d_basis01 * points[index + 1].position
                ) / tDiff;
    }

    inline InterpolationType computeCurvature(size_t index, floating_t tDiff, floating_t t) const
    {
        auto oneMinusT = 1 - t;

        //we're computing the second derivative of the computePosition function with respect to t
        //we can do this by computing the second derivatives of each of its basis functions.
        //thankfully this can easily be done analytically since they're polynomials!
        auto d2_basis00 = 60 * oneMinusT * t * (2 * t - 1);
        auto d2_basis10 = 12 * oneMinusT * t * (5 * t - 3);
        auto d2_basis20 = t * (t * (-10 * t + 18) - 9) + 1;
        auto d2_basis21 = t * (t * (10 * t - 12) + 3);
        auto d2_basis11 = 12 * oneMinusT * t * (5 * t - 2);
        auto d2_basis01 = 60 * oneMinusT * t * (1 - 2 * t);

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return (
                    d2_basis00 * points[index].position +
                    d2_basis10 * tDiff * points[index].tangent +
                    d2_basis20 * tDiff * tDiff * points[index].curvature +

                    d2_basis21 * tDiff * tDiff * points[index + 1].curvature +
                    d2_basis11 * tDiff * points[index + 1].tangent +
                    d2_basis01 * points[index + 1].position
                ) / (tDiff * tDiff);
    }

    inline InterpolationType computeWiggle(size_t index, floating_t tDiff, floating_t t) const
    {
        //we're computing the third derivative of the computePosition function with respect to t
        auto d3_basis00 = 60 * (6 * t * (1 - t) + 1);
        auto d3_basis10 = 12 * (t * (16 - 15 * t) + 3);
        auto d3_basis20 = t * (36 - 30 * t) - 9;
        auto d3_basis21 = t * (30 * t - 24) + 3;
        auto d3_basis11 = 12 * (t * (14 - 15 * t) + 2);
        auto d3_basis01 = -d3_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return (
                    d3_basis00 * points[index].position +
                    d3_basis10 * tDiff * points[index].tangent +
                    d3_basis20 * tDiff * tDiff * points[index].curvature +

                    d3_basis21 * tDiff * tDiff * points[index + 1].curvature +
                    d3_basis11 * tDiff * points[index + 1].tangent +
                    d3_basis01 * points[index + 1].position
                ) / (tDiff * tDiff * tDiff);
    }

    inline floating_t computeSegmentLength(size_t index, floating_t from, floating_t to) const
    {
        floating_t tDiff = knots[index + 1] - knots[index];
        auto segmentFunction = [this, index, tDiff](floating_t t) -> floating_t {
            auto tangent = computeTangent(index, tDiff, t);
            return tangent.length();
        };

        return tDiff * SplineLibraryCalculus::adaptiveSimpsonsIntegral(segmentFunction, from, to);
    }

private: //data
    std::vector<QuinticHermiteSplinePoint> points;
    std::vector<floating_t> knots;
};
