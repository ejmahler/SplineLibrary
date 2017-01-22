#pragma once

#include <vector>

#include "utils/spline_common.h"
#include "utils/calculus.h"

template<class InterpolationType, typename floating_t=float>
class Spline
{
public:
    Spline(std::vector<InterpolationType> originalPoints, floating_t maxT)
        :originalPoints(std::move(originalPoints)), maxT(maxT)
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
    floating_t getMaxT(void) const { return maxT; }

    const std::vector<InterpolationType> &getOriginalPoints(void) const { return originalPoints; }
    virtual bool isLooping(void) const = 0;

    //lower level functions
    virtual size_t segmentCount(void) const = 0;
    virtual size_t segmentForT(floating_t t) const = 0;
    virtual floating_t segmentT(size_t segmentIndex) const = 0;
    virtual floating_t segmentArcLength(size_t segmentIndex, floating_t a, floating_t b) const = 0;

protected:
    const floating_t maxT;

private:
    const std::vector<InterpolationType> originalPoints;
};




template<template<class, typename> class SplineCore, class InterpolationType, typename floating_t>
class SplineImpl: public Spline<InterpolationType, floating_t>
{
public:
    InterpolationType getPosition(floating_t t) const override { return common.getPosition(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t t) const override { return common.getTangent(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t t) const override { return common.getCurvature(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t t) const override { return common.getWiggle(t); }

    floating_t arcLength(floating_t a, floating_t b) const override { if(a > b) std::swap(a,b); return ArcLength::arcLength(*this,a,b); }
    floating_t totalLength(void) const override { return ArcLength::totalLength(*this); }

    bool isLooping(void) const override { return false; }

    size_t segmentCount(void) const override { return common.segmentCount(); }
    size_t segmentForT(floating_t t) const override { return common.segmentForT(t); }
    floating_t segmentT(size_t segmentIndex) const override { return common.segmentT(segmentIndex); }
    floating_t segmentArcLength(size_t segmentIndex, floating_t a, floating_t b) const override { return common.segmentLength(segmentIndex, a, b); }

protected:
    //protected constructor and destructor, so that this class can only be used as a parent class, even though it won't have any pure virtual methods
    SplineImpl(std::vector<InterpolationType> originalPoints, floating_t maxT)
        :Spline<InterpolationType, floating_t>(std::move(originalPoints), maxT)
    {}
    ~SplineImpl(void) = default;

    SplineCore<InterpolationType, floating_t> common;
};





template<template<class, typename> class SplineCore, class InterpolationType, typename floating_t>
class SplineLoopingImpl: public Spline<InterpolationType, floating_t>
{
public:
    InterpolationType getPosition(floating_t globalT) const override
    {
        floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
        return common.getPosition(wrappedT);
    }
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const override
    {
        floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
        return common.getTangent(wrappedT);
    }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const override
    {
        floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
        return common.getCurvature(wrappedT);
    }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const override
    {
        floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
        return common.getWiggle(wrappedT);
    }

    floating_t arcLength(floating_t a, floating_t b) const override
    {
        floating_t wrappedA =  SplineCommon::wrapGlobalT(a, maxT);
        floating_t wrappedB =  SplineCommon::wrapGlobalT(b, maxT);

        return ArcLength::arcLength(*this, wrappedA, wrappedB);
    }
    floating_t totalLength(void) const override { return ArcLength::totalLength(*this); }

    bool isLooping(void) const override { return true; }

    size_t segmentCount(void) const override { return common.segmentCount(); }
    size_t segmentForT(floating_t t) const override { return common.segmentForT(t); }
    floating_t segmentT(size_t segmentIndex) const override { return common.segmentT(segmentIndex); }
    floating_t segmentArcLength(size_t segmentIndex, floating_t a, floating_t b) const override { return common.segmentLength(segmentIndex, a, b); }

protected:
    //protected constructor and destructor, so that this class can only be used as a parent class, even though it won't have any pure virtual methods
    SplineLoopingImpl(std::vector<InterpolationType> originalPoints, floating_t maxT)
        :Spline<InterpolationType, floating_t>(std::move(originalPoints), maxT)
    {}
    ~SplineLoopingImpl(void) = default;

    SplineCore<InterpolationType, floating_t> common;
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
