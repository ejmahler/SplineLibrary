
#include <QtTest/QtTest>

#include "testcalculus.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);

    TestCalculus calculusTests;

    return QTest::qExec(&calculusTests, argc, argv);
}
