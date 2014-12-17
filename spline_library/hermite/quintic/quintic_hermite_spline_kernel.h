#ifndef QUINTICHERMITESPLINEKERNEL_H
#define QUINTICHERMITESPLINEKERNEL_H


class QuinticHermiteSplineKernel
{
public:
    template<class InterpolationType, typename floating_t=double>
    struct alignas(16) InterpolationData {
        //points
        InterpolationType p0;
        InterpolationType p1;

        //tangents
        InterpolationType m0;
        InterpolationType m1;

        //curvatures
        InterpolationType c0;
        InterpolationType c1;

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
                basis00 * p0 +
                basis10 * m0 +
                basis20 * c0 +
                basis21 * c1 +
                basis11 * m1 +
                basis01 * p1;
        }

        inline InterpolationType computeTangent(floating_t t) const
        {
            auto oneMinusT = 1 - t;

            //we're essentially computing the derivative of the computePosition function with respect to t
            //we can do this by computing the derivatives of each of its basis functions.
            //thankfully this can easily be done analytically since they're polynomials!
            auto basis00 = -30 * oneMinusT * oneMinusT * t * t;
            auto basis10 = -oneMinusT * oneMinusT * (3 * t - 1) * (5 * t + 1);
            auto basis20 = -0.5 * oneMinusT * oneMinusT * t * (5 * t - 2);
            auto basis21 = -0.5 * oneMinusT * t * t * (5 * t - 3);
            auto basis11 = -t * t * (3 * t - 2) * (5 * t - 6);
            auto basis01 = 30 * oneMinusT * oneMinusT * t * t;

            return (
                basis00 * p0 +
                basis10 * m0 +
                basis20 * c0 +
                basis21 * c1 +
                basis11 * m1 +
                basis01 * p1
                    ) * tDistanceInverse;
        }

        inline InterpolationType computeCurvature(floating_t t) const
        {
            auto oneMinusT = 1 - t;

            //we're essentially computing the second derivative of the computePosition function with respect to t
            //we can do this by computing the second derivatives of each of its basis functions.
            //thankfully this can easily be done analytically since they're polynomials!
            auto basis00 = 60 * oneMinusT * t * (2 * t - 1);
            auto basis10 = 12 * oneMinusT * t * (5 * t - 3);
            auto basis20 = t * (t * (-10 * t + 18) - 9) + 1;
            auto basis21 = t * (t * (10 * t - 12) + 3);
            auto basis11 = 12 * oneMinusT * t * (5 * t - 2);
            auto basis01 = -60 * oneMinusT * t * (2 * t - 1);

            return (
                basis00 * p0 +
                basis10 * m0 +
                basis20 * c0 +
                basis21 * c1 +
                basis11 * m1 +
                basis01 * p1
                    ) * (tDistanceInverse * tDistanceInverse);
        }

        inline InterpolationType computeWiggle(floating_t t) const
        {
            //we're essentially computing the third derivative of the computePosition function with respect to t
            auto basis00 = -60 * (6 * t * (t - 1) + 1);
            auto basis10 = -12 * (t * (15 * t - 16) + 3);
            auto basis20 = t * (36 - 30 * t) - 9;
            auto basis21 = t * (30 * t - 24) + 3;
            auto basis11 = -12 * (t * (15 * t - 14) + 2);
            auto basis01 = -basis00;

            return (
                basis00 * p0 +
                basis10 * m0 +
                basis20 * c0 +
                basis21 * c1 +
                basis11 * m1 +
                basis01 * p1
                    ) * (tDistanceInverse * tDistanceInverse * tDistanceInverse);
        }
    };
private:
    QuinticHermiteSplineKernel(void) = default;
};

#endif // QUINTICHERMITESPLINEKERNEL_H
