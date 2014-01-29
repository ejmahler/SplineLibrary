SplineLibrary
=============
A C++ library to collect many useful spline functions into one place.

A spline is a formula for smoothly transitioning from one data point to the next in a data set. For example, you could create a spline containing ten colors (each stored as R, G, and B values) to create a color gradient that smoothly transitions from one color to the next.

Features
-------------
* Interpolation of standard catmull-rom splines
    * Include `spline_library/cubic_hermite/cr_spline.h`, create a `CRSpline` object, and call its `getPosition` method.
    * Several more spline types. See [Spline Types](docs/SplineTypes.md) for the full list
* Looping Splines
    * To make a looping catmull-rom spline, include `spline_library/cubic_hermite/looping_cr_spline.h` and create a `LoopingCRSpline` object instead.
    * Every spline type has both looping and non-looping variants
* Compute the inverse of a spline
    * Given a data point (not necessarily on the spline, or even close to it), what T value brings the spline closest to that data point?
    * Create a SplineInverter object and call either its findClosestFast or findClosestPrecise method
* Interpolation of the first, second, and third derivatives of the spline
    * The first derivative is called the "tangent" - this is how quickly and in what direction the interpolated position is changing, per T
    * The second derivative is called the "curvature" - this is how quickly and in what direction the interpolated tangent is changing, per T
    * The third derivative is called the "wiggle" - this is how quickly and in what direction the interpolated curvature is changing, per T



Project Layout
-------------
The root of the repository is a Qt Creator project that demonstrates some uses of the library. The source for the spline code itself is in the "spline_library" directory, and the code to set up the demo is in the "demo" directory.

The demo project requires Qt 5. To build it, either run qmake with the .pro file to generate a makefile, or open the .pro file in qt Creator.

The spline_library code has no third-party dependencies, so it's safe to drop that folder directly in the source folder of your own project.

License
-------------
This code is available under the terms of the Mozilla Public License v2.0 http://mozilla.org/MPL/2.0/
