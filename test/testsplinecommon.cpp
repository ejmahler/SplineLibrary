#include "testsplinecommon.h"

#include "spline_library/vector.h"
#include "spline_library/utils/spline_common.h"

#include "common.h"

#include <vector>
#include <cmath>

#include <QtTest/QtTest>

TestSplineCommon::TestSplineCommon(QObject *parent) : QObject(parent)
{

}

void TestSplineCommon::testInnerPadding_data(void)
{
    QTest::addColumn<std::vector<Vector2>>("points");
    QTest::addColumn<int>("padding");
    QTest::addColumn<float>("alpha");
    QTest::addColumn<std::vector<float>>("expectedT");

    size_t testSize = 8;
    std::vector<Vector2> straightLine = TestDataFloat::generateStraightLineData(testSize);
    std::vector<Vector2> traingleLine = TestDataFloat::generateTriangleNumberData(testSize);

    // padding of 0
    std::vector<float> expectedEqual0 = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f };
    std::vector<float> expectedTriangle0 = { 0.0f, 0.25f, 0.75f, 1.5f, 2.5f, 3.75f, 5.25f, 7.0f };


    QTest::newRow("Equidistant, alpha=0, padding=0") << straightLine << 0 << 0.0f  << expectedEqual0;
    QTest::newRow("Equidistant, alpha=1, padding=0") << straightLine << 0 << 1.0f << expectedEqual0;
    QTest::newRow("Increasing distance, alpha=0, padding=0") << traingleLine << 0 << 0.0f << expectedEqual0;
    QTest::newRow("Increasing distance, alpha=1, padding=0") << traingleLine << 0 << 1.0f << expectedTriangle0;

    // padding of 1
    std::vector<float> expectedEqual1 = { -1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f };
    std::vector<float> expectedTriangle1 = { -0.25f, 0.0f, 0.5f, 1.25f, 2.25f, 3.5f, 5.0f, 6.75f };
    std::vector<float> expectedCentripetal1 = { -0.508553f, 0.0f, 0.719203f, 1.60004f, 2.61715f, 3.75431f, 5.0f, 6.34551f };

    QTest::newRow("Equidistant, alpha=0, padding=1") << straightLine << 1 << 0.0f  << expectedEqual1;
    QTest::newRow("Equidistant, alpha=1, padding=1") << straightLine << 1 << 1.0f << expectedEqual1;
    QTest::newRow("Increasing distance, alpha=0, padding=1") << traingleLine << 1 << 0.0f << expectedEqual1;
    QTest::newRow("Increasing distance, alpha=1, padding=1") << traingleLine << 1 << 1.0f << expectedTriangle1;

    // padding of 1, but with alpha = 0.5 instead of 1
    QTest::newRow("Equidistant, alpha=0.5, padding=1") << straightLine << 1 << 1.0f << expectedEqual1;
    QTest::newRow("Increasing distance, alpha=0.5, padding=1") << traingleLine << 1 << 0.5f << expectedCentripetal1;

    // padding of 2
    std::vector<float> expectedEqual2 = { -2.0f, -1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f };
    std::vector<float> expectedTriangle2 = { -0.75f, -0.5f, 0.0f, 0.75f, 1.75f, 3.0f, 4.5f, 6.25f };

    QTest::newRow("Equidistant, alpha=0, padding=2") << straightLine << 2 << 0.0f  << expectedEqual2;
    QTest::newRow("Equidistant, alpha=1, padding=2") << straightLine << 2 << 1.0f << expectedEqual2;
    QTest::newRow("Increasing distance, alpha=0, padding=2") << traingleLine << 2 << 0.0f << expectedEqual2;
    QTest::newRow("Increasing distance, alpha=1, padding=2") << traingleLine << 2 << 1.0f << expectedTriangle2;
}

void TestSplineCommon::testInnerPadding(void)
{
    QFETCH(std::vector<Vector2>, points);
    QFETCH(int, padding);
    QFETCH(float, alpha);
    QFETCH(std::vector<float>, expectedT);

    // make sure points and expected have the same size
    QCOMPARE(points.size(), expectedT.size());

    // compute our T values given the input parameters
    std::vector<float> actualT = SplineCommon::computeTValuesWithInnerPadding(points, alpha, padding);

    for(size_t i = 0; i < expectedT.size(); i++) {
        QCOMPARE(actualT[i], expectedT[i]);
    }
}
