#ifndef LOOPING_CR_SPLINE_H
#define LOOPING_CR_SPLINE_H

#include <vector>
#include <unordered_map>

#include "cubichermitespline.h"

class LoopingCRSpline : public CubicHermiteSpline
{
public:
    LoopingCRSpline(const std::vector<Vector3D> &points, double alpha = 0.0);
    ~LoopingCRSpline();

    Vector3D getPosition(double x) const;
    InterpolatedPT getTangent(double x) const;
    InterpolatedPTC getCurvature(double x) const;

    double getT(int index) const;
    double getMaxT(void) const;
    int getNumSegments(void) const;

    const std::vector<Vector3D> &getPoints(void) const;

    bool isLoop(void) const;

    private: //data

        double maxT;

        //original point data
        std::vector<Vector3D> points;

        //map from index to t value. it's a map and not an array so we can store negative indexes
        std::unordered_map<int,double> indexToT;

};

#endif // LOOPING_CR_SPLINE_H
