QT       += core gui opengl widgets concurrent

TARGET = SplineDemo
TEMPLATE = app

CONFIG += c++14

exists(./SplineDemo_Include.pri) {
    include(SplineDemo_Include.pri)
}

#on windows we need to manually add opengl because ???
windows: LIBS += -lopengl32

SOURCES += \
    demo/settingswidget.cpp \
    demo/settings.cpp \
    demo/mainwindow.cpp \
    demo/graphicscontroller.cpp \
    demo/benchmarker.cpp


HEADERS  += \
    demo/settingswidget.h \
    demo/settings.h \
    demo/mainwindow.h \
    demo/graphicscontroller.h \
    spline_library/spline.h \
    spline_library/utils/linearalgebra.h \
    spline_library/utils/nanoflann.hpp \
    spline_library/utils/splinesample_adaptor.h \
    demo/benchmarker.h \
    spline_library/utils/calculus.h \
    spline_library/vector.h \
    spline_library/utils/spline_common.h \
    spline_library/splines/generic_b_spline.h \
    spline_library/splines/uniform_cubic_bspline.h \
    spline_library/splines/cubic_hermite_spline.h \
    spline_library/splines/uniform_cr_spline.h \
    spline_library/splines/quintic_hermite_spline.h \
    spline_library/splines/natural_spline.h \
    spline_library/utils/arclength.h \
    spline_library/utils/splineinverter.h


FORMS    += \
    demo/settingswidget.ui \
    demo/mainwindow.ui

test {
    message(Test build)
    QT += testlib
    TARGET = UnitTests

    HEADERS += \
        test/testcalculus.h \
        test/testvector.h \
        test/testspline.h \
        test/testlinalg.h \
        test/testarclength.h \
        test/testsplinecommon.h \
        test/common.h

    SOURCES += \
        test/test_main.cpp \
        test/testcalculus.cpp \
        test/testvector.cpp \
        test/testspline.cpp \
        test/testlinalg.cpp \
        test/testarclength.cpp \
        test/testsplinecommon.cpp

} else {
    SOURCES += demo/main.cpp
}


#uncomment these to enable debug info in release mode
#QMAKE_CXXFLAGS_RELEASE += -g
#QMAKE_CFLAGS_RELEASE += -g
#QMAKE_LFLAGS_RELEASE =
