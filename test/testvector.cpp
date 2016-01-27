#include "testvector.h"

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

TestVector::TestVector(QObject *parent) : QObject(parent)
{

}

void TestVector::testConstructors(void)
{
    Vector<3, float> vec1;
    QCOMPARE(vec1[0], 0.0f);
    QCOMPARE(vec1[1], 0.0f);
    QCOMPARE(vec1[2], 0.0f);

    Vector<3, float> vec2({1.0f, 2.0f, 3.0f });
    QCOMPARE(vec2[0], 1.0f);
    QCOMPARE(vec2[1], 2.0f);
    QCOMPARE(vec2[2], 3.0f);

    Vector<3, float> vec3;
    vec3[0] = vec2[0];
    vec3[1] = vec2[1];
    vec3[2] = vec2[2];

    QVERIFY(vec2 == vec3);
    QVERIFY(vec1 != vec3);
}

void TestVector::testVectorArithmetic_data(void)
{
    QTest::addColumn<Vector3>("left");
    QTest::addColumn<Vector3>("right");
    QTest::addColumn<Vector3>("additionExpected");
    QTest::addColumn<Vector3>("subtractionExpected");

    QTest::newRow("allZero") << Vector3() << Vector3() << Vector3() << Vector3();

    QTest::newRow("sumToZero") << Vector3({1,1,1}) << Vector3({-1,-1,-1}) << Vector3() << Vector3({2,2,2});
    QTest::newRow("oneSideZero") << Vector3({1,1,1}) << Vector3() << Vector3({1,1,1}) << Vector3({1,1,1});
    QTest::newRow("nonZero") << Vector3({2,2,2}) << Vector3({1,1,1}) << Vector3({3,3,3}) << Vector3({1,1,1});
}

void TestVector::testVectorArithmetic(void)
{
    QFETCH(Vector3, left);
    QFETCH(Vector3, right);
    QFETCH(Vector3, additionExpected);
    QFETCH(Vector3, subtractionExpected);

    auto sum = left + right;
    auto reverseSum = right + left;
    QCOMPARE(sum, additionExpected);
    QCOMPARE(reverseSum, additionExpected);

    auto difference = left - right;
    auto reverseDifference = -(right - left);
    QCOMPARE(difference, subtractionExpected);
    QCOMPARE(reverseDifference, subtractionExpected);

    auto additionCopy = left;
    additionCopy += right;
    QCOMPARE(additionCopy, additionExpected);

    auto subtractionCopy = left;
    subtractionCopy -= right;
    QCOMPARE(subtractionCopy, subtractionExpected);
}

void TestVector::testScalarArithmetic_data(void)
{
    QTest::addColumn<Vector3>("left");
    QTest::addColumn<float>("right");
    QTest::addColumn<Vector3>("multiplicationExpected");
    QTest::addColumn<Vector3>("divisionExpected");

    QTest::newRow("multiplyByZero") << Vector3({1,1,1}) << 0.0f << Vector3() << Vector3();
    QTest::newRow("multiplyByOne") << Vector3({2,2,2}) << 1.0f << Vector3({2,2,2}) << Vector3({2,2,2});
    QTest::newRow("multiplyByTwo") << Vector3({2,2,2}) << 2.0f << Vector3({4,4,4}) << Vector3({1,1,1});
    QTest::newRow("multiplyByNegative") << Vector3({2,2,2}) << -2.0f << Vector3({-4,-4,-4}) << Vector3({-1,-1,-1});
}

void TestVector::testScalarArithmetic(void)
{
    QFETCH(Vector3, left);
    QFETCH(float, right);
    QFETCH(Vector3, multiplicationExpected);
    QFETCH(Vector3, divisionExpected);

    auto product = left * right;
    auto reverseProduct = right * left;
    QCOMPARE(product, multiplicationExpected);
    QCOMPARE(reverseProduct, multiplicationExpected);

    auto productCopy = left;
    productCopy *= right;
    QCOMPARE(productCopy, multiplicationExpected);

    if(right != 0) {
        auto division = left / right;
        QCOMPARE(division, divisionExpected);

        auto divisionCopy = left;
        divisionCopy /= right;
        QCOMPARE(divisionCopy, divisionExpected);
    }
}

void TestVector::testLengthOperations_data(void)
{
    QTest::addColumn<Vector3>("v");
    QTest::addColumn<float>("lengthExpected");

    QTest::newRow("zeroLength") << Vector3({0,0,0}) << 0.0f;
    QTest::newRow("oneLength") << Vector3({1,0,0}) << 1.0f;
    QTest::newRow("rootThree") << Vector3({1,1,1}) << std::sqrt(3.0f);
    QTest::newRow("pythagorean") << Vector3({3,4,12}) << 13.0f;
}

void TestVector::testLengthOperations(void)
{
    QFETCH(Vector3, v);
    QFETCH(float, lengthExpected);

    auto length = v.length();
    auto length2 = v.lengthSquared();
    QCOMPARE(length, lengthExpected);
    QCOMPARE(length2, lengthExpected * lengthExpected);

    auto normalized = v.normalized();
    auto reScaled = normalized * lengthExpected;

    //we have to compare element by element because of floating point errors
    for(size_t i = 0; i < 3; i++) {
        QCOMPARE(reScaled[i], v[i]);
    }
}

void TestVector::testSplineFunctionality_data(void)
{
    std::vector<Vector2> cubicPoints {
        Vector2({-1,-1}),
        Vector2({ 0, 0}),
        Vector2({ 1, 1}),
        Vector2({ 2, 2})
    };

    std::vector<Vector2> quinticPoints {
        Vector2({-2,-2}),
        Vector2({-1,-1}),
        Vector2({ 0, 0}),
        Vector2({ 1, 1}),
        Vector2({ 2, 2}),
        Vector2({ 3, 3})
    };

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");

    std::shared_ptr<Spline<Vector2>> uniformCubic = std::make_shared<UniformCubicBSpline<Vector2>>(cubicPoints);
    std::shared_ptr<Spline<Vector2>> genericB = std::make_shared<GenericBSpline<Vector2>>(cubicPoints, 3);
    std::shared_ptr<Spline<Vector2>> natural = std::make_shared<NaturalSpline<Vector2>>(cubicPoints,false);
    std::shared_ptr<Spline<Vector2>> uniformCR = std::make_shared<UniformCRSpline<Vector2>>(cubicPoints);
    std::shared_ptr<Spline<Vector2>> cubicHermite = std::make_shared<CubicHermiteSpline<Vector2>>(cubicPoints);
    std::shared_ptr<Spline<Vector2>> quinticHermite = std::make_shared<QuinticHermiteSpline<Vector2>>(quinticPoints);

    QTest::newRow("uniformCubicB")  << uniformCubic;
    QTest::newRow("genericB3")      << genericB;
    QTest::newRow("natural")        << natural;
    QTest::newRow("uniformCR")      << uniformCR;
    QTest::newRow("cubicHermite")   << cubicHermite;
    //QTest::newRow("quinticHermite") << quinticHermite;
}

void TestVector::testSplineFunctionality(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);

    auto positionTangent = spline->getTangent(0.5f);
    Vector2 expectedPosition({0.5f, 0.5f});
    Vector2 expectedTangent({1.0f, 1.0f});

    for(size_t i = 0; i < 2; i++) {
        QCOMPARE(positionTangent.position[i], expectedPosition[i]);
    }
    for(size_t i = 0; i < 2; i++) {
        QCOMPARE(positionTangent.tangent[i], expectedTangent[i]);
    }

    float length = spline->totalLength();
    float expectedLength = std::sqrt(2.0f);
    QCOMPARE(length, expectedLength);

    Vector2 queryPoint({0.4f, 0.0f});
    float expectedClosestT = 0.2f;

    SplineInverter<Vector2> inverter(*spline.get());
    float closestT = inverter.findClosestT(queryPoint);
    QCOMPARE(std::round(closestT * 1000)/1000, expectedClosestT);
}
