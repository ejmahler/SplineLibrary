
#include <QtTest/QtTest>

#include "testcalculus.h"
#include "testvector.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);

    TestCalculus calculusTests;
    TestVector vectorTests;

    return QTest::qExec(&calculusTests, argc, argv) | QTest::qExec(&vectorTests, argc, argv);
}
