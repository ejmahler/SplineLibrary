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
    spline_library/cubic_hermite/looping_cubic_hermite_spline.cpp \
    spline_library/cubic_hermite/looping_cr_spline.cpp \
    spline_library/cubic_hermite/cubic_hermite_spline.cpp \
    spline_library/cubic_hermite/cr_spline.cpp \
    spline_library/vector3d.cpp \
    spline_library/splineinverter.cpp \
    spline_library/spline.cpp \
    spline_library/quintic_hermite/quintic_hermite_spline.cpp \
    spline_library/quintic_hermite/quintic_cr_spline.cpp \
    spline_library/quintic_hermite/looping_quintic_hermite_spline.cpp \
    spline_library/quintic_hermite/looping_quintic_cr_spline.cpp \
    spline_library/b_spline/cubic_b_spline.cpp \
    spline_library/b_spline/looping_cubic_b_spline.cpp

HEADERS  += \
    demo/settingswidget.h \
    demo/settings.h \
    demo/mainwindow.h \
    demo/graphicscontroller.h \
    spline_library/cubic_hermite/looping_cubic_hermite_spline.h \
    spline_library/cubic_hermite/looping_cr_spline.h \
    spline_library/cubic_hermite/cubic_hermite_spline.h \
    spline_library/cubic_hermite/cr_spline.h \
    spline_library/vector3d.h \
    spline_library/splineinverter.h \
    spline_library/spline.h \
    spline_library/quintic_hermite/quintic_hermite_spline.h \
    spline_library/quintic_hermite/quintic_cr_spline.h \
    spline_library/quintic_hermite/looping_quintic_hermite_spline.h \
    spline_library/quintic_hermite/looping_quintic_cr_spline.h \
    spline_library/b_spline/cubic_b_spline.h \
    spline_library/b_spline/looping_cubic_b_spline.h

FORMS    += \
    demo/settingswidget.ui \
    demo/mainwindow.ui
