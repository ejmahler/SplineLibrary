#ifndef CUBICBSPLINEKERNEL_H
#define CUBICBSPLINEKERNEL_H

class CubicBSplineKernel
{
public:
    template<class InterpolationType, typename floating_t=double>
    struct alignas(16) InterpolationData {
        //points
        InterpolationType beforePoint, p0, p1, afterPoint;

        //t values
        floating_t t0, t1;

        inline floating_t computeLocalT(floating_t globalT) const
        {
            return globalT - t0;
        }

        inline InterpolationType computePosition(floating_t t) const
        {
            return (
                        beforePoint * ((1 - t) * (1 - t) * (1 - t)) +
                        p0 * (t * t * 3 * (t - 2) + 4) +
                        p1 * (t * (t * (-3 * t + 3) + 3) + 1) +
                        afterPoint * (t * t * t)
                    ) / 6;
        }

        inline InterpolationType computeTangent(floating_t t) const
        {
            return (
                        beforePoint * (-(1 - t) * (1 - t)) +
                        p0 * (t * (3 * t - 4)) +
                        p1 * ((3 * t + 1) * (1 - t)) +
                        afterPoint * (t * t)
                    ) / 2;
        }

        inline InterpolationType computeCurvature(floating_t t) const
        {
            return (
                        beforePoint * (1 - t) +
                        p0 * (3 * t - 2) +
                        p1 * (1 - 3 * t) +
                        afterPoint * (t)
                    );
        }

        inline InterpolationType computeWiggle(void) const
        {
            return 3 * (p0 - p1) + (afterPoint - beforePoint);
        }
    };

private:
    CubicBSplineKernel(void) = default;
};

#endif // CUBICBSPLINEKERNEL_H
