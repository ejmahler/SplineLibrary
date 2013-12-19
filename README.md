SplineLibrary
=============
A C++ library to collect many useful spline functions into one place.

A spline is a formula for smoothly transitioning from one data point to the next in a data set. For example, you could create a spline containing ten colors (each stored as R, G, and B values) to create a color gradient that smoothly transitions from one color to the next.

Features
-------------
* Interpolation of standard catmull-rom splines
    * Include `spline_source/cr_spline.h`, create a `CRSpline` object, and call its `getPosition` method.
* Looped Splines
    * To make a looped catmull-rom spline, include `spline_source/looping_cr_spline.h` and create a `LoopingCRSpline` object instead.
    * Every spline has a looped variant.
* Interpolation of chordal and centripetal catmull-rom splines
    * Set the "alpha" parameter in the constructor of CatmullRomSpline
    * 0 for a standard catmull-rom spline (default)
    * 0.5 for a centripetal catmull-rom spline
    * 1 for a chordal catmull-rom spline
    * Any number is allowed, not just these three. The demo allows you to set this freely, so feel free to experiment. For most use cases you will want either standard or centripetal.
* Interpolation of the first and second derivatives of the spline
    * The first derivative is called the "tangent" - this is how quickly and in what direction the interpolated position is changing, per T
    * The second derivative is called the "curvature" - this is how quickly and in what direction the interpolated tangent is changing, per T
* Compute the inverse of the spline
    * Given a data point that may or may not be on the spline, what T value brings the spline closest to that data point?
    * Create a SplineInverter object and call either its findClosestFast or findClosestPrecise method

Project Layout
-------------
The root of the repository is a Qt Creator project that demonstrates some uses of the library. The source for the spline code itself is in the "spline-source" directory, and the code to set up the demo is in the "demo" directory.

The demo project requires Qt, but the spline-source code has no third-party dependencies.

License
-------------
This code is available under the terms of the Mozilla Public License v2.0 http://mozilla.org/MPL/2.0/
