#pragma once

#include <cassert>

#include "../spline.h"

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

    inline size_t segmentCount(void) const
    {
        return points.size() - 1;
    }

    inline size_t segmentForT(floating_t t) const
    {
        size_t segmentIndex = SplineCommon::getIndexForT(knots, t);
        if(segmentIndex >= segmentCount())
            return segmentCount() - 1;
        else
            return segmentIndex;
    }

    inline floating_t segmentT(size_t segmentIndex) const
    {
        return knots[segmentIndex];
    }

    inline InterpolationType getPosition(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineCommon::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return computePosition(knotIndex, tDiff, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineCommon::getIndexForT(knots, globalT);
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
        size_t knotIndex = SplineCommon::getIndexForT(knots, globalT);
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
        size_t knotIndex = SplineCommon::getIndexForT(knots, globalT);
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

    inline floating_t segmentLength(size_t index, floating_t a, floating_t b) const
    {
        floating_t tDiff = knots[index + 1] - knots[index];
        auto segmentFunction = [this, index, tDiff](floating_t t) -> floating_t {
            auto tangent = computeTangent(index, tDiff, t);
            return tangent.length();
        };

        floating_t localA = (a - knots[index]) / tDiff;
        floating_t localB = (b - knots[index]) / tDiff;

        return tDiff * SplineLibraryCalculus::gaussLegendreQuadratureIntegral(segmentFunction, localA, localB);
    }

private: //methods
    inline InterpolationType computePosition(size_t index, floating_t tDiff, floating_t t) const
    {
        //this is a logical extension of the cubic hermite spline's basis functions
        //that has one basis function for t0 position, one for t1 position
        //one for t0 tangent (1st derivative of position), and one for t1 tangent
        //this adds 2 more basis functions, one for t0 curvature (2nd derivative) and t1 curvature
        //see this paper for details http://www.rose-hulman.edu/~finn/CCLI/Notes/day09.pdf
        auto basis00 = (1 - t) * (1 - t) * (1 - t) * (t * (6 * t + 3) + 1);
        auto basis10 = t * (1 - t) * (1 - t) * (1 - t) * (3 * t + 1);
        auto basis20 = floating_t(0.5) * (1 - t) * (1 - t) * (1 - t) * t * t;
        auto basis21 = floating_t(0.5) * (1 - t) * (1 - t) * t * t * t;
        auto basis11 = t * t * t * (1 - t) * (t * 3 - 4);
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
        //we're computing the derivative of the computePosition function with respect to t
        //we can do this by computing the derivatives of each of its basis functions.
        //thankfully this can easily be done analytically since they're polynomials!
        auto d_basis00 =  (-30) * (1 - t) * (1 - t) * t * t;
        auto d_basis10 = (1 - t) * (1 - t) * (1 - 3 * t) * (5 * t + 1);
        auto d_basis20 = floating_t(-0.5) * (1 - t) * (1 - t) * t * (5 * t - 2);
        auto d_basis21 = floating_t(0.5) * (1 - t) * t * t * (3 - 5 * t);
        auto d_basis11 = t * t * (2 - 3 * t) * (5 * t - 6);
        auto d_basis01 = 30 * (t - 1) * (t - 1) * t * t;

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
        //we're computing the second derivative of the computePosition function with respect to t
        //we can do this by computing the second derivatives of each of its basis functions.
        //thankfully this can easily be done analytically since they're polynomials!
        auto d2_basis00 = t * ((180 - 120 * t) * t - 60);
        auto d2_basis10 = t * ((96 - 60 * t) * t - 36);
        auto d2_basis20 = t * ((18  -10 * t) * t - 9) + 1;
        auto d2_basis21 = t * (t * (10 * t - 12) + 3);
        auto d2_basis11 = t * ((84 - 60 * t) * t - 24);
        auto d2_basis01 = -d2_basis00;

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
        auto d3_basis00 = (360 - 360*t) * t - 60;
        auto d3_basis10 = (192 - 180*t) * t - 36;
        auto d3_basis20 = (36 - 30 * t) * t - 9;
        auto d3_basis21 = (30 * t - 24) * t + 3;
        auto d3_basis11 = (168 - 180*t) * t - 24;
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

private: //data
    std::vector<QuinticHermiteSplinePoint> points;
    std::vector<floating_t> knots;
};


template<class InterpolationType, typename floating_t=float>
class QuinticHermiteSpline final : public SplineImpl<QuinticHermiteSplineCommon, InterpolationType, floating_t>
{
//constructors
public:
    QuinticHermiteSpline(const std::vector<InterpolationType> &points,
                         const std::vector<InterpolationType> &tangents,
                         const std::vector<InterpolationType> &curvatures,
                         floating_t alpha = 0.0
                         )
        :SplineImpl<QuinticHermiteSplineCommon, InterpolationType,floating_t>(points, points.size() - 1)
    {
        assert(points.size() >= 2);
        assert(points.size() == tangents.size());
        assert(points.size() == curvatures.size());

        int size = points.size();
        int padding = 0;
        int numSegments = size - 1;

        //compute the T values for each point
        auto indexToT = SplineCommon::computeTValuesWithInnerPadding(points, alpha, padding);

        //pre-arrange the data needed for interpolation
        std::vector<floating_t> knots(numSegments + 1);
        std::vector<typename QuinticHermiteSplineCommon<InterpolationType, floating_t>::QuinticHermiteSplinePoint> positionData(numSegments + 1);
        for(int i = 0; i < numSegments + 1; i++)
        {
            knots[i] = indexToT.at(i);

            positionData[i].position = points.at(i);
            positionData[i].tangent = tangents.at(i);
            positionData[i].curvature = curvatures.at(i);
        }

        common = QuinticHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
    }

    QuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha = 0.0f)
        :SplineImpl<QuinticHermiteSplineCommon, InterpolationType,floating_t>(points, points.size() - 5)
    {
        assert(points.size() >= 6);

        int size = points.size();
        int firstTangent = 1;
        int lastTangent = size - 2;
        int firstCurvature = 2;
        int numSegments = size - 5;

        //compute the T values for each point
        auto indexToT = SplineCommon::computeTValuesWithInnerPadding(points, alpha, firstCurvature);

        //compute the tangents
        std::vector<InterpolationType> tangents(size);
        for(int i = firstTangent; i <= lastTangent; i++)
        {
            floating_t tPrev = indexToT.at(i - 1);
            floating_t tCurrent = indexToT.at(i);
            floating_t tNext = indexToT.at(i + 1);

            InterpolationType pPrev = points.at(i - 1);
            InterpolationType pCurrent = points.at(i);
            InterpolationType pNext = points.at(i + 1);

            //the tangent is the standard catmull-rom spline tangent calculation
            tangents[i] =
                      pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                    + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

                 //plus a little something extra - this is derived from the pyramid contruction
                 //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
                 //yielding the standard catmull-rom formula
                    - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
        }

        //compute the curvatures
        std::vector<InterpolationType> curves(size);
        int lastCurvature = firstCurvature + numSegments;
        for(int i = firstCurvature; i <= lastCurvature; i++)
        {
            floating_t tPrev = indexToT.at(i - 1);
            floating_t tCurrent = indexToT.at(i);
            floating_t tNext = indexToT.at(i + 1);

            InterpolationType pPrev = tangents.at(i - 1);
            InterpolationType pCurrent = tangents.at(i);
            InterpolationType pNext = tangents.at(i + 1);

            //the tangent is the standard catmull-rom spline tangent calculation
            curves[i] =
                      pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                    + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

                 //plus a little something extra - this is derived from the pyramid contruction
                 //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
                 //yielding the standard catmull-rom formula
                    - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
        }


        //pre-arrange the data needed for interpolation
        int lastSegment = firstCurvature + numSegments;
        std::vector<floating_t> knots(numSegments + 1);
        std::vector<typename QuinticHermiteSplineCommon<InterpolationType, floating_t>::QuinticHermiteSplinePoint> positionData(numSegments + 1);
        for(int i = firstCurvature; i < lastSegment + 1; i++)
        {
            knots[i - firstCurvature] = indexToT.at(i);

            positionData[i - firstCurvature].position = points.at(i);
            positionData[i - firstCurvature].tangent = tangents.at(i);
            positionData[i - firstCurvature].curvature = curves.at(i);
        }
        common = QuinticHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
    }
};



template<class InterpolationType, typename floating_t=float>
class LoopingQuinticHermiteSpline final : public SplineLoopingImpl<QuinticHermiteSplineCommon, InterpolationType, floating_t>
{
//constructors
public:
    LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points,
                                const std::vector<InterpolationType> &tangents,
                                const std::vector<InterpolationType> &curvatures,
                                floating_t alpha = 0.0
                                )
        :SplineLoopingImpl<QuinticHermiteSplineCommon, InterpolationType,floating_t>(points, points.size())
    {
        assert(points.size() >= 2);
        assert(points.size() == tangents.size());
        assert(points.size() == curvatures.size());

        int size = points.size();
        int numSegments = size;

        //compute the T values for each point
        int padding = 0;
        auto indexToT = SplineCommon::computeLoopingTValues(points, alpha, padding);

        //pre-arrange the data needed for interpolation
        std::vector<floating_t> knots(numSegments + 1);
        std::vector<typename QuinticHermiteSplineCommon<InterpolationType, floating_t>::QuinticHermiteSplinePoint> positionData(numSegments + 1);
        for(int i = 0; i < numSegments + 1; i++)
        {
            knots[i] = indexToT.at(i);

            positionData[i].position = points.at(i);
            positionData[i].tangent = tangents.at(i);
            positionData[i].curvature = curvatures.at(i);
        }

        common = QuinticHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
    }

    LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha = 0.0)
        :SplineLoopingImpl<QuinticHermiteSplineCommon, InterpolationType,floating_t>(points, points.size())
    {
        assert(points.size() >= 3);

        int size = points.size();
        int numSegments = size;

        //compute the T values for each point
        int padding = 2;
        auto indexToT = SplineCommon::computeLoopingTValues(points, alpha, padding);

        //compute the tangents
        std::vector<InterpolationType> tangents(size + 1);
        for(int i = 0; i < size + 1; i++)
        {
            floating_t tPrev = indexToT.at(i - 1);
            floating_t tCurrent = indexToT.at(i);
            floating_t tNext = indexToT.at(i + 1);

            InterpolationType pPrev = points.at((i - 1 + size)%size);
            InterpolationType pCurrent = points.at((i + size)%size);
            InterpolationType pNext = points.at((i + 1)%size);

            //the tangent is the standard catmull-rom spline tangent calculation
            tangents[i] =
                      pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                    + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

                 //plus a little something extra - this is derived from the pyramid contruction
                 //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
                 //yielding the standard catmull-rom formula
                    - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
        }

        //compute the curvatures
        std::vector<InterpolationType> curves(size + 1);
        for(int i = 0; i < size + 1; i++)
        {
            floating_t tPrev = indexToT.at(i - 1);
            floating_t tCurrent = indexToT.at(i);
            floating_t tNext = indexToT.at(i + 1);

            InterpolationType pPrev = tangents.at((i - 1 + size)%size);
            InterpolationType pCurrent = tangents.at(i);
            InterpolationType pNext = tangents.at((i + 1)%size);

            //the tangent is the standard catmull-rom spline tangent calculation
            curves[i] =
                      pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                    + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

                 //plus a little something extra - this is derived from the pyramid contruction
                 //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
                 //yielding the standard catmull-rom formula
                    - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
        }


        //pre-arrange the data needed for interpolation
        std::vector<floating_t> knots(numSegments + 1);
        std::vector<typename QuinticHermiteSplineCommon<InterpolationType, floating_t>::QuinticHermiteSplinePoint> positionData(numSegments + 1);
        for(int i = 0; i < numSegments + 1; i++)
        {
            knots[i] = indexToT[i];
            positionData[i].position = points[i%size];
            positionData[i].tangent = tangents[i];
            positionData[i].curvature = curves[i];
        }

        common = QuinticHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
    }
};
