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
    spline-source/vector3d.cpp \
    spline-source/splineinverter.cpp \
    spline-source/spline.cpp \
    spline-source/quintichermitespline.cpp \
    spline-source/quintic_cr_spline.cpp \
    spline-source/looping_quintic_cr_spline.cpp \
    spline-source/cubichermitespline.cpp \
    spline-source/cr_spline.cpp \
    spline-source/looping_cr_spline.cpp

HEADERS  += \
    demo/settingswidget.h \
    demo/settings.h \
    demo/mainwindow.h \
    demo/graphicscontroller.h \
    spline-source/vector3d.h \
    spline-source/splineinverter.h \
    spline-source/spline.h \
    spline-source/quintichermitespline.h \
    spline-source/looping_quintic_cr_spline.h \
    spline-source/quintic_cr_spline.h \
    spline-source/cubichermitespline.h \
    spline-source/cr_spline.h \
    spline-source/looping_cr_spline.h

FORMS    += \
    demo/settingswidget.ui \
    demo/mainwindow.ui
