#pragma once

#include <vector>

template<class InterpolationType, typename floating_t=float>
class Spline
{

public:

    Spline(const std::vector<InterpolationType> &originalPoints)
        :originalPoints(originalPoints)
    {}

public:
    struct InterpolatedPT;

    struct InterpolatedPTC;

    struct InterpolatedPTCW;

    virtual InterpolationType getPosition(floating_t x) const = 0;
    virtual InterpolatedPT getTangent(floating_t x) const = 0;
    virtual InterpolatedPTC getCurvature(floating_t x) const = 0;
    virtual InterpolatedPTCW getWiggle(floating_t x) const = 0;

    virtual floating_t arcLength(floating_t a, floating_t b) const = 0;
    virtual floating_t totalLength(void) const = 0;

    virtual floating_t getT(int index) const = 0;
    virtual floating_t getMaxT(void) const = 0;

    const std::vector<InterpolationType> &getOriginalPoints(void) const { return originalPoints; }
    virtual bool isLooping(void) const = 0;

private:
    const std::vector<InterpolationType> originalPoints;
};

template<class InterpolationType, typename floating_t>
struct Spline<InterpolationType,floating_t>::InterpolatedPT
{
    InterpolationType position;
    InterpolationType tangent;

    InterpolatedPT(const InterpolationType &p, const InterpolationType &t)
        :position(p),tangent(t)
    {}
};

template<class InterpolationType, typename floating_t>
struct Spline<InterpolationType,floating_t>::InterpolatedPTC
{
    InterpolationType position;
    InterpolationType tangent;
    InterpolationType curvature;

    InterpolatedPTC(const InterpolationType &p, const InterpolationType &t, const InterpolationType &c)
        :position(p),tangent(t),curvature(c)
    {}
};

template<class InterpolationType, typename floating_t>
struct Spline<InterpolationType,floating_t>::InterpolatedPTCW
{
    InterpolationType position;
    InterpolationType tangent;
    InterpolationType curvature;
    InterpolationType wiggle;

    InterpolatedPTCW(const InterpolationType &p, const InterpolationType &t, const InterpolationType &c, const InterpolationType &w)
        :position(p),tangent(t),curvature(c), wiggle(w)
    {}
};
