QT       += core gui opengl widgets concurrent

TARGET = SplineDemo
TEMPLATE = app

CONFIG += c++14

#on mac we want to look for boost in homebrew folder
windows:INCLUDEPATH+= "C:\Boost\boost_1_60_0"
macx:INCLUDEPATH += /usr/local/Cellar/boost/1.59.0/include

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
    spline_library/splineinverter.h \
    spline_library/spline.h \
    spline_library/natural/looping_natural_spline.h \
    spline_library/natural/natural_spline.h \
    spline_library/hermite/cubic/cubic_hermite_spline.h \
    spline_library/hermite/cubic/looping_cubic_hermite_spline.h \
    spline_library/hermite/quintic/looping_quintic_hermite_spline.h \
    spline_library/hermite/quintic/quintic_hermite_spline.h \
    spline_library/utils/linearalgebra.h \
    spline_library/utils/nanoflann.hpp \
    spline_library/utils/splinesample_adaptor.h \
    spline_library/basis/generic_b_spline.h \
    spline_library/basis/looping_generic_b_spline.h \
    spline_library/basis/generic_b_spline_common.h \
    demo/benchmarker.h \
    spline_library/natural/natural_spline_common.h \
    spline_library/hermite/cubic/cubic_hermite_spline_common.h \
    spline_library/hermite/quintic/quintic_hermite_spline_common.h \
    spline_library/basis/uniform_cubic_bspline_common.h \
    spline_library/basis/uniform_cubic_bspline.h \
    spline_library/basis/looping_uniform_cubic_bspline.h \
    spline_library/hermite/cubic/uniform_cr_spline.h \
    spline_library/hermite/cubic/uniform_cr_spline_common.h \
    spline_library/hermite/cubic/looping_uniform_cr_spline.h \
    spline_library/utils/calculus.h \
    spline_library/vector.h \
    spline_library/utils/spline_common.h

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
        test/testlinalg.h

    SOURCES += \
        test/test_main.cpp \
        test/testcalculus.cpp \
        test/testvector.cpp \
        test/testspline.cpp \
        test/testlinalg.cpp
} else {
    SOURCES += demo/main.cpp
}


#uncomment these to enable debug info in release mode
#QMAKE_CXXFLAGS_RELEASE += -g
#QMAKE_CFLAGS_RELEASE += -g
#QMAKE_LFLAGS_RELEASE =
