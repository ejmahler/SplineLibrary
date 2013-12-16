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
    spline-source/catmullromspline.cpp

HEADERS  += \
    demo/settingswidget.h \
    demo/settings.h \
    demo/mainwindow.h \
    demo/graphicscontroller.h \
    spline-source/vector3d.h \
    spline-source/splineinverter.h \
    spline-source/spline.h \
    spline-source/catmullromspline.h

FORMS    += \
    demo/settingswidget.ui \
    demo/mainwindow.ui
