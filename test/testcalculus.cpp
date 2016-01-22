#include "testcalculus.h"

#include "spline_library/utils/calculus.h"

#include <QtTest/QtTest>

TestCalculus::TestCalculus(QObject *parent) : QObject(parent)
{

}

void TestCalculus::testAdaptiveSimpsons_data(void)
{
    QTest::addColumn<float>("from");
    QTest::addColumn<float>("to");
    QTest::addColumn<float>("expected");

    QTest::newRow("balanced") << -3.0f << 3.0f << -18.0f;
    QTest::newRow("all positive") << 2.0f << 5.0f << 113.25f;
    QTest::newRow("uneven") << -2.0f << 4.0f << 36.0f;
}

void TestCalculus::testAdaptiveSimpsons(void)
{
    QFETCH(float, from);
    QFETCH(float, to);
    QFETCH(float, expected);

    auto result = SplineLibraryCalculus::adaptiveSimpsonsIntegral([](auto x){return x*x*(x-1);}, from, to);

    QCOMPARE(result + 1, expected + 1);
}
