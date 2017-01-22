#include "testspline.h"

#include "spline_library/vector.h"
#include "spline_library/utils/spline_common.h"

#include "spline_library/splines/uniform_cubic_bspline.h"
#include "spline_library/splines/generic_b_spline.h"
#include "spline_library/splines/natural_spline.h"
#include "spline_library/splines/cubic_hermite_spline.h"
#include "spline_library/splines/uniform_cr_spline.h"
#include "spline_library/splines/quintic_hermite_spline.h"

#include "spline_library/utils/splineinverter.h"

#include "common.h"

#include <vector>
#include <memory>
#include <cmath>

#include <QtTest/QtTest>

TestSpline::TestSpline(QObject *parent) : QObject(parent)
{

}

void TestSpline::testMethods_data(void)
{
    //our data will just be points on a straight line between 0 and 100
    //this makes the total length of this line 100 * sqrt(2) so it'll be easy to verify
    std::vector<Vector2> data {
        Vector2({0,0}),
        Vector2({1,0}),
        Vector2({3,3}),
        Vector2({6,6}),
        Vector2({10,10}),
        Vector2({15,15}),
        Vector2({21,21}),
        Vector2({28,28}),
        Vector2({36,36}),
        Vector2({45,45}),
        Vector2({55,55})
    };

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<int>("padding");
    QTest::addColumn<float>("alpha");
    QTest::addColumn<size_t>("expectedSegments");

    QTest::newRow("uniformCR") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<UniformCRSpline<Vector2>>(addPadding(data,1))) << 1 << 0.0f << data.size()-1;
    QTest::newRow("catmullRom") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<CubicHermiteSpline<Vector2>>(addPadding(data,1))) << 1 << 0.0f << data.size()-1;
    QTest::newRow("catmullRomAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<CubicHermiteSpline<Vector2>>(addPadding(data,1), 0.5f)) << 1 << 0.5f << data.size()-1;

    QTest::newRow("cubicHermite") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<CubicHermiteSpline<Vector2>>(data, makeTangents(data))) << 0 << 0.0f << data.size()-1;
    QTest::newRow("cubicHermiteAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<CubicHermiteSpline<Vector2>>(data, makeTangents(data), 0.5f)) << 0 << 0.5f << data.size()-1;

    QTest::newRow("quinticCatmullRom") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<QuinticHermiteSpline<Vector2>>(addPadding(data,2))) << 2 << 0.0f << data.size()-1;
    QTest::newRow("quinticCatmullRomAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<QuinticHermiteSpline<Vector2>>(addPadding(data,2), 0.5f)) << 2 << 0.5f << data.size()-1;

    QTest::newRow("quinticHermite") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<QuinticHermiteSpline<Vector2>>(data, makeTangents(data), makeTangents(makeTangents(data)))) << 0 << 0.0f << data.size()-1;
    QTest::newRow("quinticHermiteAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<QuinticHermiteSpline<Vector2>>(data, makeTangents(data), makeTangents(makeTangents(data)), 0.5f)) << 0 << 0.5f << data.size()-1;

    QTest::newRow("natural") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(data, true))            << 0 << 0.0f << data.size()-1;
    QTest::newRow("naturalAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(data, true, 0.5f)) << 0 << 0.5f << data.size()-1;

    QTest::newRow("naturalNotAKnot") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(data, true, 0.0f, NaturalSpline<Vector2>::NotAKnot))      << 0 << 0.0f << data.size()-1;
    QTest::newRow("naturalAlphaNotAKnot") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(data, true, 0.5f, NaturalSpline<Vector2>::NotAKnot)) << 0 << 0.5f << data.size()-1;

    QTest::newRow("naturalWithoutEndpoints") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(addPadding(data,1), false))            << 1 << 0.0f << data.size()-1;
    QTest::newRow("naturalWithoutEndpointsAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(addPadding(data,1), false, 0.5f)) << 1 << 0.5f << data.size()-1;

    QTest::newRow("uniformB") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<UniformCubicBSpline<Vector2>>(addPadding(data,1)))     << 1 << 0.0f << data.size()-1;
    QTest::newRow("genericBCubic") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<GenericBSpline<Vector2>>(addPadding(data,1), 3))  << 1 << 0.0f << data.size()-1;
    QTest::newRow("genericBQuintic") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<GenericBSpline<Vector2>>(addPadding(data,2), 5)) << 2 << 0.0f << data.size()-1;
}

void TestSpline::testMethods(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(int, padding);
    QFETCH(float, alpha);
    QFETCH(size_t, expectedSegments);

    //test the methods that require no input
    float maxT = spline->getMaxT();
    QCOMPARE(maxT, float(expectedSegments));

    size_t segmentCount = spline->segmentCount();
    QCOMPARE(segmentCount, expectedSegments);

    bool looping = spline -> isLooping();
    QCOMPARE(looping, false);

    //test the "segmentT" method to make sure it returns the expected values
    auto expectedT = SplineCommon::computeTValuesWithInnerPadding(spline->getOriginalPoints(), alpha, padding);

    for(size_t i = 0; i < spline->segmentCount(); i++) {
        QCOMPARE(spline->segmentT(i), expectedT[i + padding]);
    }

    //test the "segment for T" method to make sure it returns the segment index we expect
    for(size_t i = 0; i < spline->segmentCount(); i++) {
        float beginT = expectedT[i + padding];
        QCOMPARE(spline->segmentForT(beginT), i);

        float halfwayT = (expectedT[i + padding] + expectedT[i + 1 + padding]) / 2;
        QCOMPARE(spline->segmentForT(halfwayT), i);
    }

    //make sure out-of-range T values in segmentForT return correct results
    QCOMPARE(spline->segmentForT(-10), size_t(0));
    QCOMPARE(spline->segmentForT(spline->getMaxT()), expectedSegments - 1);
    QCOMPARE(spline->segmentForT(spline->getMaxT() * 2), expectedSegments - 1);
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

    //numerically integrate the tangent
    Vector2 integratedTangent =
            gaussLegendreQuadratureIntegral(tangentFunction, spline->segmentT(0), spline->segmentT(1))
            + gaussLegendreQuadratureIntegral(tangentFunction, spline->segmentT(1), spline->segmentT(2));

    QCOMPARE(integratedTangent[0], expectedPosition[0]);
    QCOMPARE(integratedTangent[1], expectedPosition[1]);

    //numerically integrate the curvature
    Vector2 integratedCurvature =
            gaussLegendreQuadratureIntegral(curveFunction, spline->segmentT(0), spline->segmentT(1))
            + gaussLegendreQuadratureIntegral(curveFunction, spline->segmentT(1), spline->segmentT(2));

    QCOMPARE(integratedCurvature[0], expectedTangent[0]);
    QCOMPARE(integratedCurvature[1], expectedTangent[1]);

    //numerically integrate the wiggle
    //this test will fail if the spline algorithm doesn't have a continuous curvature! Spline algorithms that should not have continuous curvature
    //IE cubic hermite spline, generic b spline with degree 2) should go in the non c2 test below
    Vector2 integratedWiggle1 = gaussLegendreQuadratureIntegral(wiggleFunction, spline->segmentT(0), spline->segmentT(1));
    Vector2 integratedWiggle2 = gaussLegendreQuadratureIntegral(wiggleFunction, spline->segmentT(1), spline->segmentT(2));

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
        auto endResults = spline->getCurvature(spline->segmentT(2));
        auto midResults = spline->getCurvature(spline->segmentT(1));
        auto beginResults = spline->getCurvature(spline->segmentT(0));

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
    //numerically integrate the tangent
    Vector2 integratedTangent =
            gaussLegendreQuadratureIntegral(tangentFunction, spline->segmentT(0), spline->segmentT(1))
            + gaussLegendreQuadratureIntegral(tangentFunction, spline->segmentT(1), spline->segmentT(2));

    QCOMPARE(integratedTangent[0], expectedPosition[0]);
    QCOMPARE(integratedTangent[1], expectedPosition[1]);

    //numerically integrate the curvature
    Vector2 integratedCurvature =
            gaussLegendreQuadratureIntegral(curveFunction, spline->segmentT(0), spline->segmentT(1))
            + gaussLegendreQuadratureIntegral(curveFunction, spline->segmentT(1), spline->segmentT(2));

    QCOMPARE(integratedCurvature[0], expectedTangent[0]);
    QCOMPARE(integratedCurvature[1], expectedTangent[1]);

    //numerically integrate the wiggle. unlike the other test we're only doing the second of the two segments
    //because the code to account for the discontinuity in curvature is really awkward and non-accurate
    Vector2 integratedWiggle = gaussLegendreQuadratureIntegral(wiggleFunction, spline->segmentT(1), spline->segmentT(2));

    QCOMPARE(integratedWiggle[0], expectedCurvature[0]);
    QCOMPARE(integratedWiggle[1], expectedCurvature[1]);
}
