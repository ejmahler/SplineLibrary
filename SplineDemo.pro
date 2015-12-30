QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SplineDemo
TEMPLATE = app

CONFIG += c++14

#on mac we want to look for boost in homebrew folder
macx:INCLUDEPATH += /usr/local/Cellar/boost/1.59.0/include

SOURCES += \
    demo/settingswidget.cpp \
    demo/settings.cpp \
    demo/mainwindow.cpp \
    demo/main.cpp \
    demo/graphicscontroller.cpp

HEADERS  += \
    demo/settingswidget.h \
    demo/settings.h \
    demo/mainwindow.h \
    demo/graphicscontroller.h \
    spline_library/vector3d.h \
    spline_library/splineinverter.h \
    spline_library/spline.h \
    spline_library/basis/cubic_b_spline.h \
    spline_library/basis/looping_cubic_b_spline.h \
    spline_library/natural/looping_natural_spline.h \
    spline_library/natural/natural_spline.h \
    spline_library/hermite/cubic/cubic_hermite_spline.h \
    spline_library/hermite/cubic/looping_cubic_hermite_spline.h \
    spline_library/hermite/quintic/looping_quintic_hermite_spline.h \
    spline_library/hermite/quintic/quintic_hermite_spline.h \
    spline_library/splinelengthcalculator.h \
    spline_library/utils/linearalgebra.h \
    spline_library/utils/nanoflann.hpp \
    spline_library/utils/splinesample_adaptor.h \
    spline_library/linear/linear_spline.h \
    spline_library/linear/looping_linear_spline.h \
    spline_library/utils/spline_setup.h \
    spline_library/natural/natural_spline_kernel.h \
    spline_library/linear/linear_spline_kernel.h \
    spline_library/basis/cubic_b_spline_kernel.h \
    spline_library/hermite/cubic/cubic_hermite_spline_kernel.h \
    spline_library/hermite/quintic/quintic_hermite_spline_kernel.h \
    spline_library/basis/generic_b_spline.h

FORMS    += \
    demo/settingswidget.ui \
    demo/mainwindow.ui
