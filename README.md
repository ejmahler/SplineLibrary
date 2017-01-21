SplineLibrary
=============
A C++ library created to provide open-source reference implementations of many spline functions, e.g. Natural Splines and Catmull-Rom Splines. It's a challenge to find well-documented and FOSS implementations of these useful tools, and mathematical defintions of these things are really hard to decipher, so my hope is that this library will provide a starting point for others.

A spline is a formula for smoothly transitioning from one data point to the next in a data set. For example, you could create a spline containing ten colors (each stored as R, G, and B values) to create a color gradient that smoothly transitions from one color to the next.

Features
-------------
* Interpolation of catmull-rom splines
    * Include `spline_library/hermite/cubic/uniform_cr_spline.h`, create a `UniformCRSpline` object, and call its `getPosition` method.
    * Several more spline types. See [Spline Types](docs/SplineTypes.md) for the full list
* Looping Splines
	* Also called "Periodic" or "Cyclic": These splines form a loop, in that the ending connects with the beginning
	* Calling getPosition(t) with an out-of-range T value will _wrap around" to the other end of the spline
    * To make a looping catmull-rom spline, include `spline_library/hermite/cubic/looping_uniform_cr_spline.h` and create a `LoopingUniformCRSpline` object instead.
    * Every spline type has both looping and non-looping variants
* Compute the arc length of a spline
	* Call a spline's totalLength() method to find the arc length of the entire spline
    * Call a spline's arcLength(a,b) method to find the arc length between two arbitrary T values
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

[Spline Utilities](docs/SplineUtilities.md) - Documentation of some utility classes for splines.

Project Layout
-------------
The root of the repository is a Qt Creator project that demonstrates some uses of the library. The source for the spline code itself is in the "spline_library" directory, and the code to set up the demo is in the "demo" directory.

Usage
-------------
Drop the spline_library directory in the root source folder of your project. It's header-only, so from here all you need to do is import it from your own code.

`spline_library/spline_inverter.h` and `spline_library/arclength.h` depend on Boost's Math module. If you don't want to install Boost, you can safely avoid including these two files - nothing else includes it, and nothing else relies on Boost.

Both the demo and the spline_library code require a fully compliant C++14 compiler.

Demo
-------------
Follow these steps to run the demo (Assuming you already have Qt 5.5+ installed and working):

1. Install Boost. On Linux, this is in most package managers. On Mac, it can be installed via homebrew. Otherwise, visit http://www.boost.org/
2. Create a file called `SplineDemo_Include.pri` in the root of the project
3. In this file, paste the following, where `/path/to/boost` contains Boost's include files.

    `INCLUDEPATH += "/path/to/boost"`
     
    On windows, this might be:  
    `INCLUDEPATH+= "C:\Boost\boost_1_60_0"`
     
    On Mac, this might be:
    INCLUDEPATH += /usr/local/Cellar/boost/1.59.0/include
4. Run qmake on `SplineDemo.pro` to generate a makefile, then build the makefile, and run the compiled executable
5. OR, open `SplineDemo.pro` in Qt Creator and press play

License
-------------
This code is available under the [Simplified BSD License](http://opensource.org/licenses/BSD-2-Clause)

This project includes the [nanoflann](https://github.com/jlblancoc/nanoflann) library for fast nearest-neighbor queries, which is also available under the [Simplified BSD License](http://opensource.org/licenses/BSD-2-Clause)

To-Do
-------------
* Implement "composite splines", ie a spline that is made by combining two or more splines. Current ideas include a "sum spline" where the output for a given T is a sum of all the child spline results at T, and a "concatenation spline" formed simply by starting one spline where the previous leaves off.
* More spline types as I discover them
* Find an actual mathematical definition for the quintic catmull-rom spline. The quintic cubic hermite spline is well-defined, but I basically guessed on how to automatically compute the tangents and curvatures based on the input points for the catmull-rom equivalent.
