SplineLibrary
=============
A C++ library to collect many useful spline functions into one place.

Features
-------------
* Interpolation of standard catmull-rom splines
    * Simply create a CatmullRomSpline object and call its getPosition method.
* Looped Splines
    * To make a looped spline, add a final point to the data set that is exactly equal to the first
* Interpolation of chordal and centripetal catmull-rom splines
    * Set the "alpha" parameter in the constructor of CatmullRomSpline
    * 0.5 for a centripetal catmull-rom spline (default)
    * 1 for a chordal catmull-rom spline
    * 0 for a standard catmull-rom spline
* Interpolation of the first and second derivatives of the spline
    * For lack of a better name, the first derivative is called the "velocity" - this is how quickly and in what direction the interpolated position is changing, per T
    * For lack of a better name, the second derivative is called the "acceleration" - this is how quickly and in what direction the interpolated velocity is changing, per T
* Compute the inverse of the spline
    * Given a data point that may or may not be on the spline, what T value brings the spline closest to that data point?
    * Create a SplineInverter object and call either its findClosestFast or findClosestPrecise method

Project Layout
-------------
The root of the repository is a Qt Creator project that demonstrates some uses of the library. The source for the spline code itself is in the "spline-source" directory, and the code to set up the demo is in the "demo" directory.

License
-------------
This code is available under the terms of the Mozilla Public License v2.0 http://mozilla.org/MPL/2.0/
