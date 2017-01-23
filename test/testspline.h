#pragma once

#include <array>
#include <QObject>
#include "spline_library/vector.h"

class TestSpline : public QObject
{
    Q_OBJECT
public:
    explicit TestSpline(QObject *parent = 0);

signals:

private slots:

    //test each spline's basic methods like maxT, segmentCount
    void testMethods_data(void);
    void testMethods(void);

    //use numeric integration to verify the first, second, and third derivatives
    void testDerivatives_data(void);
    void testDerivatives(void);

    //Verify that the 'segment arc length' method computes the correct result
    void testSegmentArcLength_data(void);
    void testSegmentArcLength(void);
};
