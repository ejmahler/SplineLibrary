#include "testcalculus.h"

#include "spline_library/utils/calculus.h"

#include <QtTest/QtTest>

TestCalculus::TestCalculus(QObject *parent) : QObject(parent)
{

}

void TestCalculus::testGaussLegendre_data(void)
{
    QTest::addColumn<float>("from");
    QTest::addColumn<float>("to");
    QTest::addColumn<float>("expected");

    QTest::newRow("balanced") << -3.0f << 3.0f << -18.0f;
    QTest::newRow("no transformation") << -1.0f << 1.0f << -2.0f / 3.0f;
    QTest::newRow("all positive") << 2.0f << 5.0f << 113.25f;
    QTest::newRow("uneven") << -2.0f << 4.0f << 36.0f;
}

void TestCalculus::testGaussLegendre(void)
{
    QFETCH(float, from);
    QFETCH(float, to);
    QFETCH(float, expected);

    auto result = SplineLibraryCalculus::gaussLegendreQuadratureIntegral<float>([](auto x){return x*x*(x-1);}, from, to);

    QCOMPARE(result, expected);
}
