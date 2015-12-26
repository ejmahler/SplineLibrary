SplineLibrary
=============
A C++ library created to provide open-source reference implementations of many spline functions, e.g. Natural Splines and Catmull-Rom Splines. It's a challenge to find well-documented and FOSS implementations of these useful tools, and mathematical defintions of these things are really hard to decipher, so my hope is that this library will provide a starting point for others.

A spline is a formula for smoothly transitioning from one data point to the next in a data set. For example, you could create a spline containing ten colors (each stored as R, G, and B values) to create a color gradient that smoothly transitions from one color to the next.

Features
-------------
* Interpolation of catmull-rom splines
    * Include `spline_library/hermite/cubic/cubic_hermite_spline.h`, create a `CubicHermiteSpline` object, and call its `getPosition` method.
    * Several more spline types. See [Spline Types](docs/SplineTypes.md) for the full list
* Looping Splines
    * To make a looping catmull-rom spline, include `spline_library/hermite/cubic/looping_cubic_hermite_spline.h` and create a `LoopingCubicHermiteSpline` object instead.
    * Every spline type has both looping and non-looping variants
* Compute the inverse of a spline
    * Given a data point (not necessarily on the spline, or even close to it), what T value brings the spline closest to that data point?
    * Create a SplineInverter object and call its findClosestT method
* Computation of the first, second, and third derivatives of the spline
    * The first derivative is called the "tangent" - this is how quickly and in what direction the interpolated position is changing, per T
    * The second derivative is called the "curvature" - this is how quickly and in what direction the interpolated tangent is changing, per T
    * The third derivative is called the "wiggle" - this is how quickly and in what direction the interpolated curvature is changing, per T
    

Documentation
-------------
[Glossary](docs/Glossary.md) - Glossary of important terms for understanding splines.

[Spline class API](docs/SplineAPI.md) - API documentation of the `Spline` base class.

[Spline Types](docs/SplineTypes.md) - Complete list of all supported spline formulas

[Spline Utilities](docs/SplineUtilities.md) - Documentation of some utility classes for splines, including `SplineInverter` and `SplineLengthCalculator`.

Project Layout
-------------
The root of the repository is a Qt Creator project that demonstrates some uses of the library. The source for the spline code itself is in the "spline_library" directory, and the code to set up the demo is in the "demo" directory.

Requirements
-------------

Both the demo and the spline_library code require a fully compliant C++14 compiler.

The demo project requires Qt 5.5. To build it, either run qmake with the .pro file to generate a makefile, or open the .pro file in qt Creator.

When actually using splines in your own project, drop the spline_library directory in the root source folder. It's header-only, so from here all you need to do is import it from your own code.

`spline_library/spline_inverter.h` depends on Boost's Math module. If you don't want to install Boost and you don't need the spline inverter's functionality, you can safely avoid including `spline_library/spline_inverter.h` - nothing else includes it, and nothing else relies on Boost.

There are no other third-party dependencies within the spline_library directory.

License
-------------
This code is available under the [Simplified BSD License](http://opensource.org/licenses/BSD-2-Clause)

This project includes the [nanoflann](https://github.com/jlblancoc/nanoflann) library for fast nearest-neighbor queries, which is also available under the [Simplified BSD License](http://opensource.org/licenses/BSD-2-Clause)

To-Do
-------------
* Implement a Generic B-Sline type that supports arbitrary powers, in addition to the cubic-only version. The generic one would almost certainly have worse performance and numerical stability than the cubic-only version, so it's worth keeping both.
* Add support for arbitrary T value differences in the Cubic B-Spline
* Implement "composite splines", ie a spline that is made by combining two or more splines. Current ideas include a "sum spline" where the output for a given T is a sum of all the child spline results at T, and a "concatenation spline" formed simply by starting one spline where the previous leaves off.
