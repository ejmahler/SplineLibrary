#ifndef GENERIC_B_SPLINE_COMMON
#define GENERIC_B_SPLINE_COMMON

#include <vector>

template<class InterpolationType, typename floating_t>
class GenericBSplineCommon
{
public:
    inline GenericBSplineCommon(void) = default;
    inline GenericBSplineCommon(std::vector<InterpolationType> positions, std::vector<floating_t> knots, int splineDegree)
        :positions(std::move(positions)), knots(std::move(knots)), splineDegree(splineDegree)
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);
        return computeDeboor(startIndex + 1, splineDegree, globalT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 2)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 2),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 3)
                    );
    }

//methods
private:
    InterpolationType computeDeboor(size_t knotIndex, int degree, float globalT) const;
    InterpolationType computeDeboorDerivative(size_t knotIndex, int degree, float globalT, int derivativeLevel) const;

private:
    std::vector<InterpolationType> positions;
    std::vector<floating_t> knots;
    int splineDegree;
};

template<class InterpolationType, typename floating_t>
InterpolationType GenericBSplineCommon<InterpolationType,floating_t>::computeDeboor(size_t knotIndex, int degree, float globalT) const
{
    if(degree == 0)
    {
        return positions[knotIndex];
    }
    else
    {
        floating_t alpha = (globalT - knots[knotIndex - 1]) / (knots[knotIndex + splineDegree - degree] - knots[knotIndex - 1]);

        InterpolationType leftRecursive = computeDeboor(knotIndex - 1, degree - 1, globalT);
        InterpolationType rightRecursive = computeDeboor(knotIndex, degree - 1, globalT);

        InterpolationType blended = leftRecursive * (1 - alpha) + rightRecursive * alpha;

        return blended;
    }
}

template<class InterpolationType, typename floating_t>
InterpolationType GenericBSplineCommon<InterpolationType,floating_t>::computeDeboorDerivative(size_t knotIndex, int degree, float globalT, int derivativeLevel) const
{
    if(degree == 0)
    {
        return InterpolationType();
    }
    else
    {
        float multiplier = degree / (knots[knotIndex + splineDegree - degree] - knots[knotIndex - 1]);

        if(derivativeLevel <= 1)
        {
            //once we reach this point we don't want to compute the derivative anymore
            return multiplier *
                    (computeDeboor(knotIndex, degree - 1, globalT)
                   - computeDeboor(knotIndex - 1, degree - 1, globalT)
                     );
        }
        else
        {
            //recursively call the derivative function to compute a higher derivative
            return multiplier *
                    (computeDeboorDerivative(knotIndex, degree - 1, globalT, derivativeLevel - 1)
                   - computeDeboorDerivative(knotIndex - 1, degree - 1, globalT, derivativeLevel - 1)
                     );
        }
    }
}

#endif // GENERIC_B_SPLINE_COMMON

