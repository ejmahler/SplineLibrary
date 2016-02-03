#include "testspline.h"

#include "spline_library/vector.h"

#include "spline_library/basis/uniform_cubic_bspline.h"
#include "spline_library/basis/generic_b_spline.h"
#include "spline_library/natural/natural_spline.h"
#include "spline_library/hermite/cubic/cubic_hermite_spline.h"
#include "spline_library/hermite/cubic/uniform_cr_spline.h"
#include "spline_library/hermite/quintic/quintic_hermite_spline.h"

#include "spline_library/splineinverter.h"

#include <vector>

#include <memory>

#include <QtTest/QtTest>

Q_DECLARE_METATYPE(Vector2)
Q_DECLARE_METATYPE(Vector3)
Q_DECLARE_METATYPE(std::shared_ptr<Spline<Vector2>>)

TestSpline::TestSpline(QObject *parent) : QObject(parent)
{

}

void TestSpline::testDerivatives_data(void)
{
    std::vector<Vector2> cubicPoints {
        Vector2({-4,-1}),
        Vector2({ 0, 1}),
        Vector2({ 1, 3}),
        Vector2({ 6, -4}),
        Vector2({ 5, 0}),
    };

    std::vector<Vector2> quinticPoints {
        Vector2({-2,-2}),
        Vector2({-4,-1}),
        Vector2({ 0, 1}),
        Vector2({ 2, 3}),
        Vector2({ 1, 1}),
        Vector2({ 2, 1}),
        Vector2({ 3, 2}),
    };

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<Vector2>("expectedPosition");
    QTest::addColumn<Vector2>("expectedTangent");
    QTest::addColumn<Vector2>("expectedCurvature");

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {
        auto endResults = spline->getCurvature(spline->getMaxT());
        auto beginResults = spline->getCurvature(0);

        QTest::newRow(name)
                << spline
                << endResults.position - beginResults.position
                << endResults.tangent - beginResults.tangent
                << endResults.curvature - beginResults.curvature;
    };

    rowFunction("uniformCubicB", std::make_shared<UniformCubicBSpline<Vector2>>(cubicPoints));
    rowFunction("genericB3", std::make_shared<GenericBSpline<Vector2>>(cubicPoints, 3));
    rowFunction("natural", std::make_shared<NaturalSpline<Vector2>>(cubicPoints,false));
    rowFunction("naturalAlpha1", std::make_shared<NaturalSpline<Vector2>>(cubicPoints, false, 1.0f));
    rowFunction("quinticHermite", std::make_shared<QuinticHermiteSpline<Vector2>>(quinticPoints));
    rowFunction("quinticHermiteAlpha1", std::make_shared<QuinticHermiteSpline<Vector2>>(quinticPoints, 1.0f));
}

void TestSpline::testDerivatives(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(Vector2, expectedPosition);
    QFETCH(Vector2, expectedTangent);
    QFETCH(Vector2, expectedCurvature);

    auto tangentFunction = [=](float t) {
        return spline->getTangent(t).tangent;
    };
    auto curveFunction = [=](float t) {
        return spline->getCurvature(t).curvature;
    };
    auto wiggleFunction = [=](float t) {
        return spline->getWiggle(t).wiggle;
    };

    //find the "first" t in the spline
    size_t firstT = 0;
    while(spline->getT(firstT) < 0)
        firstT++;

    //numerically integrate the tangent
    Vector2 integratedTangent =
            gaussLegendreQuadratureIntegral(tangentFunction, spline->getT(firstT), spline->getT(firstT + 1))
            + gaussLegendreQuadratureIntegral(tangentFunction, spline->getT(firstT + 1), spline->getT(firstT + 2));

    QCOMPARE(integratedTangent[0], expectedPosition[0]);
    QCOMPARE(integratedTangent[1], expectedPosition[1]);

    //numerically integrate the curvature
    Vector2 integratedCurvature =
            gaussLegendreQuadratureIntegral(curveFunction, spline->getT(firstT), spline->getT(firstT + 1))
            + gaussLegendreQuadratureIntegral(curveFunction, spline->getT(firstT + 1), spline->getT(firstT + 2));

    QCOMPARE(integratedCurvature[0], expectedTangent[0]);
    QCOMPARE(integratedCurvature[1], expectedTangent[1]);

    //numerically integrate the wiggle
    //this test will fail if the spline algorithm doesn't have a continuous curvature! Spline algorithms that should not have continuous curvature
    //IE cubic hermite spline, generic b spline with degree 2) should go in the non c2 test below
    Vector2 integratedWiggle1 = gaussLegendreQuadratureIntegral(wiggleFunction, spline->getT(firstT), spline->getT(firstT + 1));
    Vector2 integratedWiggle2 = gaussLegendreQuadratureIntegral(wiggleFunction, spline->getT(firstT + 1), spline->getT(firstT + 2));

    Vector2 integratedWiggle = integratedWiggle1 + integratedWiggle2;

    QCOMPARE(integratedWiggle[0], expectedCurvature[0]);
    QCOMPARE(integratedWiggle[1], expectedCurvature[1]);
}



void TestSpline::testDerivativesNonC2_data(void)
{
    std::vector<Vector2> cubicPoints {
        Vector2({-4,-1}),
        Vector2({ 0, 1}),
        Vector2({ 1, 3}),
        Vector2({ 6, -4}),
        Vector2({ 5, 0}),
    };

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<Vector2>("expectedPosition");
    QTest::addColumn<Vector2>("expectedTangent");
    QTest::addColumn<Vector2>("expectedCurvature");

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {
        auto endResults = spline->getCurvature(spline->getMaxT());
        auto midResults = spline->getCurvature(spline->getT(2));
        auto beginResults = spline->getCurvature(0);

        QTest::newRow(name)
                << spline
                << endResults.position - beginResults.position
                << endResults.tangent - beginResults.tangent
                << endResults.curvature - midResults.curvature;
    };
    rowFunction("uniformCR", std::make_shared<UniformCRSpline<Vector2>>(cubicPoints));
    rowFunction("cubicHermite", std::make_shared<CubicHermiteSpline<Vector2>>(cubicPoints));
    rowFunction("cubicHermiteAlpha1", std::make_shared<CubicHermiteSpline<Vector2>>(cubicPoints, 1.0f));
}

void TestSpline::testDerivativesNonC2(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(Vector2, expectedPosition);
    QFETCH(Vector2, expectedTangent);
    QFETCH(Vector2, expectedCurvature);

    auto tangentFunction = [=](float t) {
        return spline->getTangent(t).tangent;
    };
    auto curveFunction = [=](float t) {
        return spline->getCurvature(t).curvature;
    };
    auto wiggleFunction = [=](float t) {
        return spline->getWiggle(t).wiggle;
    };

    //find the "first" t in the spline
    size_t firstT = 0;
    while(spline->getT(firstT) < 0)
        firstT++;

    //numerically integrate the tangent
    Vector2 integratedTangent =
            gaussLegendreQuadratureIntegral(tangentFunction, spline->getT(firstT), spline->getT(firstT + 1))
            + gaussLegendreQuadratureIntegral(tangentFunction, spline->getT(firstT + 1), spline->getT(firstT + 2));

    QCOMPARE(integratedTangent[0], expectedPosition[0]);
    QCOMPARE(integratedTangent[1], expectedPosition[1]);

    //numerically integrate the curvature
    Vector2 integratedCurvature =
            gaussLegendreQuadratureIntegral(curveFunction, spline->getT(firstT), spline->getT(firstT + 1))
            + gaussLegendreQuadratureIntegral(curveFunction, spline->getT(firstT + 1), spline->getT(firstT + 2));

    QCOMPARE(integratedCurvature[0], expectedTangent[0]);
    QCOMPARE(integratedCurvature[1], expectedTangent[1]);

    //numerically integrate the wiggle. unlike the other test we're only doing the second of the two segments
    //because the code to account for the discontinuity in curvature is really awkward and non-accurate
    Vector2 integratedWiggle = gaussLegendreQuadratureIntegral(wiggleFunction, spline->getT(firstT + 1), spline->getT(firstT + 2));

    QCOMPARE(integratedWiggle[0], expectedCurvature[0]);
    QCOMPARE(integratedWiggle[1], expectedCurvature[1]);
}

