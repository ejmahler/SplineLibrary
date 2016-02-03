
#include <QtTest/QtTest>

#include "testcalculus.h"
#include "testvector.h"
#include "testspline.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);

    TestCalculus calculusTests;
    TestVector vectorTests;
    TestSpline splineTests;

    return QTest::qExec(&calculusTests, argc, argv) | QTest::qExec(&vectorTests, argc, argv) | QTest::qExec(&splineTests, argc, argv);
}
