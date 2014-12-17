#ifndef LINEARSPLINEKERNEL_H
#define LINEARSPLINEKERNEL_H


class LinearSplineKernel
{
public:
    template<class InterpolationType, typename floating_t=double>
    struct alignas(16) InterpolationData
    {
        //points
        InterpolationType p0, p1;

        //t values
        floating_t t0, t1;

        //reciprocal of distance in T between p0 and p1
        floating_t tDistanceInverse;

        inline floating_t computeLocalT(floating_t globalT) const
        {
            return (globalT - t0) * tDistanceInverse;
        }

        inline InterpolationType computePosition(floating_t t)
        {
            return p0 * (1 - t) + p1 * t;
        }

        inline InterpolationType computeTangent(void)
        {
            return (p1 - p0) * tDistanceInverse;
        }
    };
private:
    LinearSplineKernel(void) = default;
};

#endif // LINEARSPLINEKERNEL_H
