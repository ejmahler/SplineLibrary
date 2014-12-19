#ifndef CUBICHERMITESPLINEKERNEL_H
#define CUBICHERMITESPLINEKERNEL_H


class CubicHermiteSplineKernel
{
public:
    template<class InterpolationType, typename floating_t>
    struct alignas(16) InterpolationData {
        //points
        InterpolationType p0;
        InterpolationType p1;

        //tangents
        InterpolationType m0;
        InterpolationType m1;

        //t values
        floating_t t0;
        floating_t t1;

        //reciprocal of distance in T between p0 and p1
        floating_t tDistanceInverse;

        inline floating_t computeLocalT(floating_t globalT) const
        {
            return (globalT - t0) * tDistanceInverse;
        }

        inline InterpolationType computePosition(floating_t t) const
        {
            auto oneMinusT = 1 - t;

            auto basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
            auto basis10 = t * oneMinusT * oneMinusT;

            auto basis11 = t * t * -oneMinusT;
            auto basis01 = t * t * (3 - 2*t);

            return
                    basis00 * p0 +
                    basis10 * m0 +

                    basis11 * m1 +
                    basis01 * p1;
        }

        inline InterpolationType computeTangent(floating_t t) const
        {
            auto oneMinusT = 1 - t;

            auto d_basis00 = -6 * t * oneMinusT;
            auto d_basis10 = -(3*t - 1) * oneMinusT;

            auto d_basis11 = t * (3 * t - 2);
            auto d_basis01 = -d_basis00;

            //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
            //intuitively it would just be the derivative of the position function and nothing else
            //if you know why please let me know
            return (
                    d_basis00 * p0 +
                    d_basis10 * m0 +

                    d_basis11 * m1 +
                    d_basis01 * p1
                    ) * tDistanceInverse;
        }

        inline InterpolationType computeCurvature(floating_t t) const
        {
            auto d2_basis00 = 6 * (2 * t - 1);
            auto d2_basis10 = 2 * (3 * t - 2);

            auto d2_basis11 = 2 * (3 * t - 1);
            auto d2_basis01 = -d2_basis00;

            //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
            //intuitively it would just be the 2nd derivative of the position function and nothing else
            //if you know why please let me know
            return (
                    d2_basis00 * p0 +
                    d2_basis10 * m0 +

                    d2_basis11 * m1 +
                    d2_basis01 * p1
                    ) * (tDistanceInverse * tDistanceInverse);
        }

        inline InterpolationType computeWiggle(void) const
        {
            //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
            //intuitively it would just be the 2nd derivative of the position function and nothing else
            //if you know why please let me know
            return (
                        12 * (p0 - p1) + 6 * (m0 + m1)
                    ) * (tDistanceInverse * tDistanceInverse * tDistanceInverse);
        }
    };

private:
    CubicHermiteSplineKernel(void) = default;
};

#endif // CUBICHERMITESPLINEKERNEL_H
