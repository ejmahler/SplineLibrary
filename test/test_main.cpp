
#include <QtTest/QtTest>

#include "testcalculus.h"
#include "testvector.h"
#include "testspline.h"
#include "testlinalg.h"
#include "testarclength.h"
#include "testsplinecommon.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);

    TestCalculus calculusTests;
    TestVector vectorTests;
    TestSpline splineTests;
    TestLinAlg algebraTests;
    TestArcLength lengthTests;
    TestSplineCommon commonTests;

    return QTest::qExec(&calculusTests, argc, argv)
            | QTest::qExec(&vectorTests, argc, argv)
            | QTest::qExec(&splineTests, argc, argv)
            | QTest::qExec(&algebraTests, argc, argv)
            | QTest::qExec(&lengthTests, argc, argv)
            | QTest::qExec(&commonTests, argc, argv);
}
