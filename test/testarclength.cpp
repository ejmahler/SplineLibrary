#include "testarclength.h"

#include "common.h"

#include "spline_library/basis/uniform_cubic_bspline.h"
#include "spline_library/basis/generic_b_spline.h"
#include "spline_library/natural/natural_spline.h"
#include "spline_library/hermite/cubic/cubic_hermite_spline.h"
#include "spline_library/hermite/cubic/uniform_cr_spline.h"
#include "spline_library/hermite/quintic/quintic_hermite_spline.h"

#include <QtTest/QtTest>

TestArcLength::TestArcLength(QObject *parent) : QObject(parent)
{

}

void TestArcLength::testArcLengthTotalLength_data(void)
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

void TestArcLength::testArcLengthTotalLength(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);

    float arc = spline->arcLength(0, spline->getMaxT());
    float total = spline->totalLength();

    QCOMPARE(arc, total);
}


void TestArcLength::testKnownArcLength_data(void)
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

void TestArcLength::testKnownArcLength(void)
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

void TestArcLength::testSolve_data(void)
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

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {

        //add a row for just part of the spline. we want to make sure a and b fall partway through a segment
        //so we'll explicitly get segment boundaries via spline->getT and lerp between them
        float partialA = lerp(spline->segmentT(1), spline->segmentT(2), 0.75f);
        float partialB = lerp(spline->segmentT(spline->segmentCount() - 3), spline->segmentT(spline->segmentCount() - 2), 0.25f);
        std::string partialName = QString("%1 (Partial)").arg(name).toStdString();
        QTest::newRow(partialName.data()) << spline << partialA << partialB;

        //add a row where a and b are in the same segment, since this is a special case
        size_t testIndex = 3;
        float sameSegmentA = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.2f);
        float sameSegmentB = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.6f);
        std::string sameSegmentName = QString("%1 (Same)").arg(name).toStdString();
        QTest::newRow(sameSegmentName.data()) << spline << sameSegmentA << sameSegmentB;
    };

    rowFunction("uniformCR", std::make_shared<UniformCRSpline<Vector2>>(addPadding(data,1)));
    rowFunction("cubicHermiteAlpha", std::make_shared<CubicHermiteSpline<Vector2>>(addPadding(data,1), 0.5f));
}

void TestArcLength::testSolve(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(float, a);
    QFETCH(float, b);

    float arcLength = ArcLength::arcLength(*spline.get(), a, b);

    float calculatedB = ArcLength::solveLength(*spline.get(), a, arcLength);
    QCOMPARE(calculatedB, b);

    //verify that if the desiredLength is longer than the spline, maxT is returned
    float totalLength = ArcLength::totalLength(*spline.get());
    float calculatedOverLength = ArcLength::solveLength(*spline.get(), a, totalLength);
    QCOMPARE(calculatedOverLength, spline->getMaxT());
}
