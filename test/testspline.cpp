#include "testspline.h"

#include "spline_library/vector.h"
#include "spline_library/utils/spline_common.h"
#include "spline_library/utils/calculus.h"

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
    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<int>("padding");
    QTest::addColumn<float>("alpha");
    QTest::addColumn<size_t>("expectedSegments");

    //our data will just be points on a straight line
    //this makes the total arc length of this line (data[size-1] - data[0]).length() so it'll be easy to verify
    auto data = TestDataFloat::generateTriangleNumberData(10);

    QTest::newRow("uniformCR") <<       TestDataFloat::createUniformCR(data)        << 1 << 0.0f << data.size()-1;
    QTest::newRow("catmullRom") <<      TestDataFloat::createCatmullRom(data, 0.0f) << 1 << 0.0f << data.size()-1;
    QTest::newRow("catmullRomAlpha") << TestDataFloat::createCatmullRom(data, 0.5f) << 1 << 0.5f << data.size()-1;

    QTest::newRow("cubicHermite") <<        TestDataFloat::createCubicHermite(data, 0.0f) << 0 << 0.0f << data.size()-1;
    QTest::newRow("cubicHermiteAlpha") <<   TestDataFloat::createCubicHermite(data, 0.5f) << 0 << 0.5f << data.size()-1;

    QTest::newRow("quinticCatmullRom") <<       TestDataFloat::createQuinticCatmullRom(data, 0.0f) << 2 << 0.0f << data.size()-1;
    QTest::newRow("quinticCatmullRomAlpha") <<  TestDataFloat::createQuinticCatmullRom(data, 0.5f) << 2 << 0.5f << data.size()-1;

    QTest::newRow("quinticHermite") <<      TestDataFloat::createQuinticHermite(data, 0.0f) << 0 << 0.0f << data.size()-1;
    QTest::newRow("quinticHermiteAlpha") << TestDataFloat::createQuinticHermite(data, 0.5f) << 0 << 0.5f << data.size()-1;

    QTest::newRow("natural") <<         TestDataFloat::createNatural(data, true, 0.0f) << 0 << 0.0f << data.size()-1;
    QTest::newRow("naturalAlpha") <<    TestDataFloat::createNatural(data, true, 0.5f) << 0 << 0.5f << data.size()-1;

    QTest::newRow("naturalNotAKnot") <<         TestDataFloat::createNotAKnot(data, true, 0.0f) << 0 << 0.0f << data.size()-1;
    QTest::newRow("naturalAlphaNotAKnot") <<    TestDataFloat::createNotAKnot(data, true, 0.5f) << 0 << 0.5f << data.size()-1;

    QTest::newRow("naturalWithoutEndpoints") <<         TestDataFloat::createNatural(data, false, 0.0f) << 1 << 0.0f << data.size()-1;
    QTest::newRow("naturalWithoutEndpointsAlpha") <<    TestDataFloat::createNatural(data, false, 0.5f) << 1 << 0.5f << data.size()-1;

    QTest::newRow("uniformB") << TestDataFloat::createUniformBSpline(data) << 1 << 0.0f << data.size()-1;

    QTest::newRow("genericBCubic") <<   TestDataFloat::createGenericBSpline(data, 3) << 1 << 0.0f << data.size()-1;
    QTest::newRow("genericBQuintic") << TestDataFloat::createGenericBSpline(data, 5) << 2 << 0.0f << data.size()-1;
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
    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");

    auto data = TestDataFloat::generateRandomData(8);

    QTest::newRow("uniformCubicB") <<       TestDataFloat::createUniformBSpline(data);
    QTest::newRow("genericB3") <<           TestDataFloat::createGenericBSpline(data,3);
    QTest::newRow("natural") <<             TestDataFloat::createNatural(data, true, 0.0f);
    QTest::newRow("naturalAlpha") <<        TestDataFloat::createNatural(data, true, 0.5f);
    QTest::newRow("quinticHermite") <<      TestDataFloat::createQuinticHermite(data, 0.0f);
    QTest::newRow("quinticHermiteAlpha1")<< TestDataFloat::createQuinticHermite(data, 0.5f);
    QTest::newRow("UniformCR") <<           TestDataFloat::createUniformCR(data);
    QTest::newRow("cubicHermite") <<        TestDataFloat::createCubicHermite(data, 0.0f);
    QTest::newRow("cubicHermiteAlpha") <<   TestDataFloat::createCubicHermite(data, 0.5f);
}

void TestSpline::testDerivatives(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);

    auto tangentFunction = [=](float t) {
        return spline->getTangent(t).tangent;
    };
    auto curveFunction = [=](float t) {
        return spline->getCurvature(t).curvature;
    };
    auto wiggleFunction = [=](float t) {
        return spline->getWiggle(t).wiggle;
    };

    //test the tangent within each segment
    for(size_t i = 0; i < spline->segmentCount(); i++)
    {
        Vector2 expectedPosition = spline->getPosition(spline->segmentT(i+1)) - spline->getPosition(spline->segmentT(i));
        Vector2 integratedTangent = SplineLibraryCalculus::gaussLegendreQuadratureIntegral<Vector2>(tangentFunction, spline->segmentT(i), spline->segmentT(i+1));

        compareFloatsLenient(integratedTangent[0], expectedPosition[0], 0.001f);
        compareFloatsLenient(integratedTangent[1], expectedPosition[1], 0.001f);
    }

    //test the curvature within each segment
    for(size_t i = 0; i < spline->segmentCount(); i++)
    {
        Vector2 expectedTangent = tangentFunction(spline->segmentT(i+1)) - tangentFunction(spline->segmentT(i));
        Vector2 integratedCurvature = SplineLibraryCalculus::gaussLegendreQuadratureIntegral<Vector2>(curveFunction, spline->segmentT(i), spline->segmentT(i+1));

        compareFloatsLenient(integratedCurvature[0], expectedTangent[0], 0.001f);
        compareFloatsLenient(integratedCurvature[1], expectedTangent[1], 0.001f);
    }

    //test the wiggle within each segment
    //note the "-.0001f" for T values - this makes sure the T value we're testing falls inside the segment we want
    //this is critical for spline types where the curvature is not continuous, because the curvature is expected to jump at the segment boundary
    //and we don't want that jump reflected in our expected result
    for(size_t i = 0; i < spline->segmentCount(); i++)
    {
        Vector2 expectedCurvature = curveFunction(spline->segmentT(i+1)-.0001f) - curveFunction(spline->segmentT(i));
        Vector2 integratedWiggle = SplineLibraryCalculus::gaussLegendreQuadratureIntegral<Vector2>(wiggleFunction, spline->segmentT(i), spline->segmentT(i+1)-.0001f);

        compareFloatsLenient(integratedWiggle[0], expectedCurvature[0], 0.001f);
        compareFloatsLenient(integratedWiggle[1], expectedCurvature[1], 0.001f);
    }
}





void TestSpline::testSegmentArcLength_data(void)
{
    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");

    auto data = TestDataFloat::generateTriangleNumberData(10);

    QTest::newRow("uniformCubicB") <<       TestDataFloat::createUniformBSpline(data);
    QTest::newRow("genericB3") <<           TestDataFloat::createGenericBSpline(data,3);
    QTest::newRow("natural") <<             TestDataFloat::createNatural(data, true, 0.0f);
    QTest::newRow("naturalAlpha") <<        TestDataFloat::createNatural(data, true, 0.5f);
    QTest::newRow("quinticHermite") <<      TestDataFloat::createQuinticHermite(data, 0.0f);
    QTest::newRow("quinticHermiteAlpha1")<< TestDataFloat::createQuinticHermite(data, 0.5f);
    QTest::newRow("UniformCR") <<           TestDataFloat::createUniformCR(data);
    QTest::newRow("cubicHermite") <<        TestDataFloat::createCubicHermite(data, 0.0f);
    QTest::newRow("cubicHermiteAlpha") <<   TestDataFloat::createCubicHermite(data, 0.5f);
}

void TestSpline::testSegmentArcLength(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);

    auto arcLengthDerivative = [=](float t) {
        return spline->getTangent(t).tangent.length();
    };

    auto arcLength2ndDerivative = [=](float t) {
        auto interpolationResult = spline->getCurvature(t);
        return Vector2::dotProduct(interpolationResult.tangent.normalized(), interpolationResult.curvature);
    };

    for(size_t i = 0; i < spline->segmentCount(); i++)
    {
        //test the whole segment
        float segmentBeginT = spline->segmentT(i);
        float segmentEndT = spline->segmentT(i+1);

        float expectedLength = (spline->getPosition(segmentBeginT) - spline->getPosition(segmentEndT)).length();
        float actualLength = spline->segmentArcLength(i, segmentBeginT, segmentEndT);

        QCOMPARE(actualLength, expectedLength);

        //test just part of the segment
        float partialBeginT = lerp(segmentBeginT, segmentEndT, 0.25f);
        float partialEndT = lerp(segmentBeginT, segmentEndT, 0.75f);

        float expectedPartialLength = (spline->getPosition(partialBeginT) - spline->getPosition(partialEndT)).length();
        float actualPartialLength = spline->segmentArcLength(i, partialBeginT, partialEndT);

        QCOMPARE(actualPartialLength, expectedPartialLength);

        //verify the 1st derivative
        float integratedTangentLength = SplineLibraryCalculus::gaussLegendreQuadratureIntegral<float>(arcLengthDerivative, segmentBeginT, segmentEndT);
        QCOMPARE(integratedTangentLength, expectedLength);

        //verify the 2nd derivative
        float integrated2ndDerivative = SplineLibraryCalculus::gaussLegendreQuadratureIntegral<float>(arcLength2ndDerivative, segmentBeginT, segmentEndT);
        float expected2ndDerivativeResult = arcLengthDerivative(segmentEndT) - arcLengthDerivative(segmentBeginT);
        compareFloatsLenient(integrated2ndDerivative, expected2ndDerivativeResult, 0.0001f);
    }
}
