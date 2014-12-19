#ifndef NATURALSPLINEKERNEL_H
#define NATURALSPLINEKERNEL_H

class NaturalSplineKernel
{
public:
    template<class InterpolationType, typename floating_t>
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

        inline InterpolationType computePosition(floating_t t)
        {
            return a + t * (b + t * (c + t * d));
        }

        inline InterpolationType computeTangent(floating_t t)
        {
            //compute the derivative of the position function
            return b + t * (2 * c + (3 * t) * d);
        }

        inline InterpolationType computeCurvature(floating_t t)
        {
            //compute the 2nd derivative of the position function
            return 2 * c + (6 * t) * d;
        }

        inline InterpolationType computeWiggle(void)
        {
            //compute the 3rd derivative of the position function
            return 6 * d;
        }
    };
private:
    NaturalSplineKernel(void) = default;
};

#endif // NATURALSPLINEKERNEL_H
