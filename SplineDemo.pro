QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SplineDemo
TEMPLATE = app

CONFIG += c++11

SOURCES += \
    demo/settingswidget.cpp \
    demo/settings.cpp \
    demo/mainwindow.cpp \
    demo/main.cpp \
    demo/graphicscontroller.cpp \
    spline_library/vector3d.cpp \
    spline_library/splineinverter.cpp \
    spline_library/spline.cpp \
    spline_library/utils/linearsolver.cpp \
    spline_library/utils/t_calculator.cpp \
    spline_library/basis/cubic_b_spline.cpp \
    spline_library/basis/looping_cubic_b_spline.cpp \
    spline_library/natural/natural_spline.cpp \
    spline_library/natural/looping_natural_spline.cpp \
    spline_library/hermite/cubic/cr_spline.cpp \
    spline_library/hermite/cubic/cubic_hermite_spline.cpp \
    spline_library/hermite/cubic/looping_cr_spline.cpp \
    spline_library/hermite/cubic/looping_cubic_hermite_spline.cpp \
    spline_library/hermite/quintic/looping_quintic_cr_spline.cpp \
    spline_library/hermite/quintic/looping_quintic_hermite_spline.cpp \
    spline_library/hermite/quintic/quintic_cr_spline.cpp \
    spline_library/hermite/quintic/quintic_hermite_spline.cpp \
    spline_library/basis/generic_b_spline.cpp

HEADERS  += \
    demo/settingswidget.h \
    demo/settings.h \
    demo/mainwindow.h \
    demo/graphicscontroller.h \
    spline_library/vector3d.h \
    spline_library/splineinverter.h \
    spline_library/spline.h \
    spline_library/utils/linearsolver.h \
    spline_library/utils/t_calculator.h \
    spline_library/basis/cubic_b_spline.h \
    spline_library/basis/looping_cubic_b_spline.h \
    spline_library/natural/looping_natural_spline.h \
    spline_library/natural/natural_spline.h \
    spline_library/hermite/cubic/cr_spline.h \
    spline_library/hermite/cubic/cubic_hermite_spline.h \
    spline_library/hermite/cubic/looping_cr_spline.h \
    spline_library/hermite/cubic/looping_cubic_hermite_spline.h \
    spline_library/hermite/quintic/looping_quintic_cr_spline.h \
    spline_library/hermite/quintic/looping_quintic_hermite_spline.h \
    spline_library/hermite/quintic/quintic_cr_spline.h \
    spline_library/hermite/quintic/quintic_hermite_spline.h \
    spline_library/basis/generic_b_spline.h

FORMS    += \
    demo/settingswidget.ui \
    demo/mainwindow.ui
