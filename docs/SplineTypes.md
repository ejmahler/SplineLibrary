Spline Types
=============
This page gives a breakdown of each spline type, how to use each one, and the advantages/disadvantages of each type.

Tl;dr
-------------
If you need the spline to pass through the input points, start with a [Catmull-Rom Spline](#catmull-rom-spline). If you don't, start with a [Cubic B-Spline](#cubic-b-spline).

Simple Types
------------- 
If you're not sure which one to use, start with these three.

### Natural Spline
The Natural Spline computes the curvature for each point, based on the position of every point, then interpolates the spline based on the list of points and the corresponding list of curvatures.

To use, import the appropriate header:
`#include "spline_library/natural/natural_spline.h"`

Create a Natural Spline by passing a std::vector<Vector3D> to the constructor, containing a list of points to interpolate through:
`std::shared_ptr<Spline> mySpline = std::make_shared<NaturalSpline<QVector2D>>(myPointList);`

##### Advantages
* Curvature is continuous [(?)](Glossary.md#continuous-curvature)

##### Disadvantages
* No local control [(?)](Glossary.md#local-control)

### Catmull-Rom Spline
A Catmull-Rom Spline computes the tangent for each point from the positions of the two closest points, then interpolates based on both the position and the tangent.

To use, import the appropriate header:
`#include "spline_library/hermite/cubic/cubic_hermite_spline.h"`

Create a catmull-rom spline by passing a std::vector<Vector3D> to the constructor, containing a list of points to interpolate through:
`std::shared_ptr<Spline> mySpline = std::make_shared<CubicHermiteSpline<QVector2D>>(myPointList);`

##### Advantages
* Local control [(?)](Glossary.md#local-control)

##### Disadvantages
* Curvature isn't continuous [(?)](Glossary.md#continuous-curvature). For some use cases this isn't a problem, so I wouldn't worry about it unless you know you need it to be continuous.
* There **must** be a nonzero distance between each adjacent set of points
* Non-looping variation requires an "extra" point on either end of the data set which will not be interpolated

### Cubic B-Spline
The B-Spline (Basis Spline) is very similar in concept to the Bezier Curve, and the cubic B-Spline is a specific type of B-Spline.

It is possible to create B-Splines with arbitrary powers (as opposed to enforcing cubic) but enforcing cubic allows for much simpler formulas and better performance.

To use, import the appropriate header:
`#include "spline_library/basis/cubic_b_spline.h"`

Create a Cubic B-Spline by passing a std::vector<Vector3D> to the constructor, containing a list of control points:
`std::shared_ptr<Spline> mySpline = std::make_shared<CubicBSpline<QVector2D>>(myPointList);`

##### Advantages
* Local control [(?)](Glossary.md#local-control)
* Curvature is continuous [(?)](Glossary.md#continuous-curvature)

##### Disadvantages
* The interpolated line does not necessarily pass through the specified points
* Non-looping variation requires an "extra" point on either end of the data set which will not be interpolated

Advanced Types
-------------
If one of the simple types above doesn't meet your needs, the following types are also available. These spline types are more powerful, but also more difficult to use correctly, and/or carry important caveats that may not be immediately obvious.

### Centripetal Catmull-Rom Spline
The Centripetal CR Spline is a variation of the Catmull-Rom Spline formula. Instead of spacing each point exactly one T apart, the distance in T between any two points will be proportional to the square root of distance between the two points. Thus, points that are very far apart will be further apart in T than points that are close together.

To use it, provide a value for the optional `alpha` parameter in the `CubicHermiteSpline` constructor. A value of 0.5 will produce a centripetal Catmull-Rom Spline, while a value of 0.0 (default) will revert to the standard formula. Other values are allowed too - a value of 1.0 will result in a "chordal" variation, and the formula will work with any number, negative or positive. Values other than 0.0 or 0.5 should be very rare, however.

It has been proven mathematically that the centripetal variation avoids certain types of self-intersections, cusps, and overshoots, producing a more aesthetically pleasing spline.

##### Advantages (compared to CRSpline)
* Proven to avoid self-intersections and overshoots when there are large variations in distance between adjacent points.

##### Disadvantages (compared to CRSpline)
* Modifies T values of points - points that are close together will have a smaller T distance and vice versa. This may be a problem if the points are keyframes for an animation, for example, or any other data series where the T values have some external meaning

### Cubic Hermite Spline
The Cubic Hermite Spline takes a list of points, and a corresponding list of tangents for each point. The Catmull-Rom Spline is a special type of the Cubic Hermite Spline which automatically computes the tangents, rather than expecting the user to supply them.

An example use case for this spline type is for physical simulation time series data, where spline->getPosition(t) returns the object's position at time T. If you know the object's velocity in addition to its position, you can make the interpolation more accurate by providing that velocity as the tangent.

To use, import the appropriate header:
`#include "spline_library/hermite/cubic/cubic_hermite_spline.h"`

Create a Cubic Hermite Spline by passing two equal-length std::vector<Vector3D> to the constructor, one containing a list of points to interpolate through, and the other containing the corresponding tangent for each point:
`std::shared_ptr<Spline> mySpline = std::make_shared<CubicHermiteSpline<QVector2D>>(myPointList, myTangentList);`

##### Advantages
* Local control [(?)](Glossary.md#local-control)
* Easily control the tangent at each point

##### Disadvantages
* Curvature isn't continuous [(?)](Glossary.md#continuous-curvature). For some use cases this isn't a problem, so I wouldn't worry about it unless you know you need it to be continuous.
* You cannot create a spline where there is zero distance between two adjacent points
* Cannot be used if you don't know the desired tangent for each point

### Quintic Hermite Spline
The Quintic Hermite Spline takes a list of points, a corresponding list of tangents for each point, and a corresponding list of curvatures for each point.

An example use case for this spline type is for physical simulation time series data, where spline->getPosition(t) returns the object's position at time T. If you know the object's acceleration in addition to its velocity and position, you can make the interpolation more accurate by providing that acceleration for the curvature.

To use, import the appropriate header:
`#include "spline_library/hermite/quintic/quintic_hermite_spline.h"`

Create a Quintic Hermite Spline by passing three equal-length std::vector<Vector3D> to the constructor, one containing a list of points to interpolate through, another containing the corresponding tangent for each point, and a third containing the corresponding curvature for each point:
`std::shared_ptr<Spline> mySpline = std::make_shared<QuinticHermiteSpline<QVector2D>>(myPointList, myTangentList, myCurvatureList);`

##### Advantages
* Local control [(?)](Glossary.md#local-control)
* Easily control the tangent and curvature at each point
* Curvature is continuous [(?)](Glossary.md#continuous-curvature)

##### Disadvantages
* You cannot create a spline where there is zero distance between two adjacent points
* Cannot be used if you don't know the desired tangent and curvature for each point
* More computationally intensive than the cubic version
* More "wiggly" than the cubic version. This sounds vague, but it's actually quantifiable: For the cubic version, the derivative of curvature is constant, but for the quintic version, the derivative of curvature is a quadractic function.
