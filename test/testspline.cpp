#include "testspline.h"

#include "spline_library/vector.h"
#include "spline_library/utils/spline_common.h"

#include "spline_library/basis/uniform_cubic_bspline.h"
#include "spline_library/basis/generic_b_spline.h"
#include "spline_library/natural/natural_spline.h"
#include "spline_library/hermite/cubic/cubic_hermite_spline.h"
#include "spline_library/hermite/cubic/uniform_cr_spline.h"
#include "spline_library/hermite/quintic/quintic_hermite_spline.h"

#include "spline_library/splineinverter.h"

#include <vector>
#include <memory>
#include <cmath>

#include <QtTest/QtTest>

Q_DECLARE_METATYPE(Vector2)
Q_DECLARE_METATYPE(Vector3)
Q_DECLARE_METATYPE(std::shared_ptr<Spline<Vector2>>)


//perform a linear interpolation between a and b
template<class InterpolationType, typename floating_t>
InterpolationType lerp(InterpolationType a, InterpolationType b, floating_t t) {
    return a * (floating_t(1) - t) + b * t;
}

//we need to pad out the ends of the data differently depending on spline type
//this way all of the splines will have the same arc length, so it'll be easier to test
template<class T>
std::vector<T> addPadding(std::vector<T> list, size_t paddingSize) {
    list.reserve(list.size() + paddingSize * 2);
    for(size_t i = 0; i < paddingSize; i++) {
        list.insert(list.begin(), list[0] - (list[1] - list[0]));
    }
    for(size_t i = 0; i < paddingSize; i++) {
        list.push_back(list[list.size() - 1] + (list[list.size() - 1] - list[list.size() - 2]));
    }
    return list;
};

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
    QTest::newRow("cubicHermite") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<CubicHermiteSpline<Vector2>>(addPadding(data,1))) << 1 << 0.0f << data.size()-1;
    QTest::newRow("cubicHermiteAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<CubicHermiteSpline<Vector2>>(addPadding(data,1), 0.5f)) << 1 << 0.5f << data.size()-1;

    QTest::newRow("quinticHermite") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<QuinticHermiteSpline<Vector2>>(addPadding(data,2))) << 2 << 0.0f << data.size()-1;
    QTest::newRow("quinticHermiteAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<QuinticHermiteSpline<Vector2>>(addPadding(data,2), 0.5f)) << 2 << 0.5f << data.size()-1;

    QTest::newRow("natural") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(data, true))            << 0 << 0.0f << data.size()-1;
    QTest::newRow("naturalAlpha") << std::static_pointer_cast<Spline<Vector2>>(std::make_shared<NaturalSpline<Vector2>>(data, true, 0.5f)) << 0 << 0.5f << data.size()-1;

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



void TestSpline::testArcLengthTotalLength_data(void)
{
    std::vector<Vector2> data {
        Vector2({100,100}),
        Vector2({400,100}),
        Vector2({500,400}),
        Vector2({300,600}),
        Vector2({300,300}),
        Vector2({150,200}),
        Vector2({100,400})
    };

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {
        QTest::newRow(name) << spline;
    };

    rowFunction("uniformCR", std::make_shared<UniformCRSpline<Vector2>>(data));
    rowFunction("cubicHermite", std::make_shared<CubicHermiteSpline<Vector2>>(data));
    rowFunction("cubicHermiteAlpha", std::make_shared<CubicHermiteSpline<Vector2>>(data, 0.5f));

    rowFunction("quinticHermite", std::make_shared<QuinticHermiteSpline<Vector2>>(data));
    rowFunction("quinticHermiteAlpha", std::make_shared<QuinticHermiteSpline<Vector2>>(data, 0.5f));

    rowFunction("natural", std::make_shared<NaturalSpline<Vector2>>(data, true));
    rowFunction("naturalAlph1", std::make_shared<NaturalSpline<Vector2>>(data, true, 0.5f));

    rowFunction("uniformB", std::make_shared<UniformCubicBSpline<Vector2>>(data));
    rowFunction("genericBCubic", std::make_shared<GenericBSpline<Vector2>>(data, 3));
    rowFunction("genericBQuintic", std::make_shared<GenericBSpline<Vector2>>(data, 5));
}

void TestSpline::testArcLengthTotalLength(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);

    float arc = spline->arcLength(0, spline->getMaxT());
    float total = spline->totalLength();

    QCOMPARE(arc, total);
}


void TestSpline::testKnownArcLength_data(void)
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
    QTest::addColumn<float>("a");
    QTest::addColumn<float>("b");
    QTest::addColumn<float>("expectedLength");
    QTest::addColumn<bool>("sameSegment");
    QTest::addColumn<size_t>("segmentIndex");

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {

        //add a test row for the whole spline
        std::string allName = QString("%1 (All)").arg(name).toStdString();
        QTest::newRow(allName.data()) << spline << 0.0f << spline->getMaxT() << (data[data.size() - 1] - data[0]).length() << false << size_t(0);

        //add a row for just part of the spline. we want to make sure a and b fall partway through a segment
        //so we'll explicitly get segment boundaries via spline->getT and lerp between them
        float partialA = lerp(spline->segmentT(2), spline->segmentT(3), 0.75f);
        float partialB = lerp(spline->segmentT(spline->segmentCount() - 3), spline->segmentT(spline->segmentCount() - 2), 0.25f);
        std::string partialName = QString("%1 (Partial)").arg(name).toStdString();
        QTest::newRow(partialName.data()) << spline << partialA << partialB << (spline->getPosition(partialA) - spline->getPosition(partialB)).length() << false << size_t(0);

        //add a row where a and b are in the same segment, since most splines treat this as a special case
        size_t testIndex = 3;
        float sameSegmentA = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.2f);
        float sameSegmentB = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.6f);
        std::string sameSegmentName = QString("%1 (Same)").arg(name).toStdString();
        QTest::newRow(sameSegmentName.data()) << spline << sameSegmentA << sameSegmentB << (spline->getPosition(sameSegmentA) - spline->getPosition(sameSegmentB)).length() << true << testIndex;
    };

    rowFunction("uniformCR", std::make_shared<UniformCRSpline<Vector2>>(addPadding(data,1)));
    rowFunction("cubicHermite", std::make_shared<CubicHermiteSpline<Vector2>>(addPadding(data,1)));
    rowFunction("cubicHermiteAlpha", std::make_shared<CubicHermiteSpline<Vector2>>(addPadding(data,1), 0.5f));

    rowFunction("quinticHermite", std::make_shared<QuinticHermiteSpline<Vector2>>(addPadding(data,2)));
    rowFunction("quinticHermiteAlpha", std::make_shared<QuinticHermiteSpline<Vector2>>(addPadding(data,2), 0.5f));

    rowFunction("natural", std::make_shared<NaturalSpline<Vector2>>(data, true));
    rowFunction("naturalAlpha", std::make_shared<NaturalSpline<Vector2>>(data, true, 0.5f));

    rowFunction("uniformB", std::make_shared<UniformCubicBSpline<Vector2>>(addPadding(data,1)));
    rowFunction("genericBCubic", std::make_shared<GenericBSpline<Vector2>>(addPadding(data,1), 3));
    rowFunction("genericBQuintic", std::make_shared<GenericBSpline<Vector2>>(addPadding(data,2), 5));
}

void TestSpline::testKnownArcLength(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(float, a);
    QFETCH(float, b);
    QFETCH(float, expectedLength);
    QFETCH(bool, sameSegment);

    float arcResult = spline->arcLength(a, b);

    auto compareFloatsLenient = [] (auto actual, auto expected) {
        auto error = std::abs(actual - expected) / expected;
        if(error > 0.01) {
            std::string errorMessage = QString("Compared floats were different. Actual: %1, Expected: %2").arg(QString::number(actual), QString::number(expected)).toStdString();
            QFAIL(errorMessage.data());
        }
    };

    //qt's fuzzy compare is a little too strict here. this is an inherently imprecise operation (especially given the use of the spline inverter)
    //so we need to allow for small deviations
    compareFloatsLenient(arcResult, expectedLength);

    //if a and b are in the same segment, we also need to test the "segmentArcLength" method for the tested segment
    if(sameSegment) {
        QFETCH(size_t, segmentIndex);

        float segmentBeginT = spline->segmentT(segmentIndex);
        float segmentEndT = spline->segmentT(segmentIndex + 1);

        float aPercent = (a - segmentBeginT) / (segmentEndT - segmentBeginT);
        float bPercent = (b - segmentBeginT) / (segmentEndT - segmentBeginT);

        float segmentArc = spline->segmentArcLength(segmentIndex, aPercent, bPercent);

        //no need for a lenient comparison here, because the results should be much closer to identical to the "arcLength" result
        QCOMPARE(arcResult, segmentArc);
    }
}
