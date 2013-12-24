SplineLibrary
=============
A C++ library to collect many useful spline functions into one place.

A spline is a formula for smoothly transitioning from one data point to the next in a data set. For example, you could create a spline containing ten colors (each stored as R, G, and B values) to create a color gradient that smoothly transitions from one color to the next.

Features
-------------
* Interpolation of standard catmull-rom splines
    * Include `spline_library/cubic_hermite/cr_spline.h`, create a `CRSpline` object, and call its `getPosition` method.
* Compute the inverse of a spline
    * Given a data point (not necessarily on the spline, or even close to it), what T value brings the spline closest to that data point?
    * Create a SplineInverter object and call either its findClosestFast or findClosestPrecise method
* Interpolation of the first and second derivatives of the spline
    * The first derivative is called the "tangent" - this is how quickly and in what direction the interpolated position is changing, per T
    * The second derivative is called the "curvature" - this is how quickly and in what direction the interpolated tangent is changing, per T
* Looping Splines
    * To make a looping catmull-rom spline, include `spline_library/cubic_hermite/looping_cr_spline.h` and create a `LoopingCRSpline` object instead.
    * Every spline type has both looping and non-looping variants
* Centripetal catmull-rom splines
    * Set the "alpha" parameter in the constructor of CRSpline
    * 0 for a standard catmull-rom spline (default)
    * 0.5 for a centripetal catmull-rom spline
    * 1 for a chordal catmull-rom spline
    * Any number is allowed, not just these three. The demo allows you to set this freely, so feel free to experiment. For most use cases you will want either standard or centripetal.
* More spline types
    * Raw cubic hermite splines
        * Import `spline_library/cubic_hermite/cubic_hermite_spline.h` and create a `CubicHermiteSpline` object
        * In addition to a list of points, provide a corresponding list of tangents
        * `CRSpline` is a subclass of `CubicHermiteSpline` which automatically computes the tangent list based on the points
    * Quintic catmull-rom splines
        * Quintic version of the catmull-rom spline
        * Import `spline_library/quintic_hermite/quintic_cr_spline.h` and create a `QuinticCRSpline` object
    * Raw quintic hermite splines
        * Import `spline_library/quintic_hermite/quintic_hermite_spline.h` and create a `QuinticHermiteSpline` object
        * In addition to a list of points, provide a corresponding list of tangents and a corresponding list of curvatures
        * `QuinticCRSpline` is a subclass of `QuinticHermiteSpline` which automatically computes the tangent list and curvature list based on the points



Project Layout
-------------
The root of the repository is a Qt Creator project that demonstrates some uses of the library. The source for the spline code itself is in the "spline_library" directory, and the code to set up the demo is in the "demo" directory.

The demo project requires Qt 5, but the spline_library code has no third-party dependencies.

License
-------------
This code is available under the terms of the Mozilla Public License v2.0 http://mozilla.org/MPL/2.0/
