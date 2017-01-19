#pragma once

#include <QObject>

class TestArcLength : public QObject
{
    Q_OBJECT
public:
    explicit TestArcLength(QObject *parent = 0);

private slots:
    //Verify arcLength(0,maxT) equals totalLength() for all non-looping splines
    void testArcLengthTotalLength_data(void);
    void testArcLengthTotalLength(void);

    //For some known arc length values, verify that each spline gives the correct result
    void testKnownArcLength_data(void);
    void testKnownArcLength(void);

    //verify that the "solve arc length" method works as expected
    void testSolve_data(void);
    void testSolve(void);

    //verify that the "partition" method works as expected
    void testPartition_data(void);
    void testPartition(void);

    //verify that the "partitionN" method works as expected
    void testPartitionN_data(void);
    void testPartitionN(void);
};
