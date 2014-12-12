#ifndef NATURALSPLINEKERNEL_H
#define NATURALSPLINEKERNEL_H

namespace NaturalSplineKernel
{
    template<class InterpolationType, typename floating_t=double>
    struct alignas(16) InterpolationData
    {
        //for natural splines, the "point data" is coefficcients for a single cubic function
        InterpolationType a, b, c, d;

        //t values
        floating_t t0;
        floating_t t1;

        //reciprocal of distance in T between t0 and t1
        floating_t tDistanceInverse;

        inline floating_t computeLocalT(floating_t globalT) const
        {
            return globalT - t0;
        }
    };

    template<class InterpolationType, typename floating_t>
    inline InterpolationType computePosition(floating_t t, const InterpolationData<InterpolationType,floating_t> &segment)
    {
        return segment.a + t * (segment.b + t * (segment.c + t * segment.d));
    }

    template<class InterpolationType, typename floating_t>
    inline InterpolationType computeTangent(floating_t t, const InterpolationData<InterpolationType,floating_t> &segment)
    {
        //compute the derivative of the position function
        return segment.b + t * (2 * segment.c + (3 * t) * segment.d);
    }

    template<class InterpolationType, typename floating_t>
    inline InterpolationType computeCurvature(floating_t t, const InterpolationData<InterpolationType,floating_t> &segment)
    {
        //compute the 2nd derivative of the position function
        return 2 * segment.c + (6 * t) * segment.d;
    }

    template<class InterpolationType, typename floating_t>
    inline InterpolationType computeWiggle(const InterpolationData<InterpolationType,floating_t> &segment)
    {
        //compute the 3rd derivative of the position function
        return 6 * segment.d;
    }

};

#endif // NATURALSPLINEKERNEL_H
