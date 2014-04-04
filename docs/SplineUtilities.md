This file outlines some utility classes for interacting with splines.

Spline Inverter
=============
The Spline Inverter object is used to answer the question "Given a data point (called the Query Point), what t value brings the spline closest to that point?" This can be used for query points which are on the spline itself, but it can also be used for points that are near the spline, or even for query points that are far away.

An example use case is an AI system for a race track. The spline could present the optimal line through a curve in the track, and the AI system could use its current position as a query point to the Spline Inverter to get the closest point on the spline. This will tell the AI where it should be, what direction it should be going via the tangent, how the optimal line is turning via the curvature, etc.

This must be computed numerically - there is no reasonable analytic solution like there is for derivatives.

To create a Spline Inverter, import `spline_library/splineinverter.h`, and create a new `SplineInverter` object by passing a `std::shared_ptr<Spline>` to the constructor. `SplineInverter inverter(mySpline);`.

In the SplineInverter constructor, it takes "samples" of the spline at regular intervals. By default it takes 10 samples per T, but this can be changed via a constructor parameter. When given a query point, it first finds the closest sample to the query point, then uses that sample location as the starting point for a refining algorithm.

### double findClosestT(const Vector3D &queryPoint) const
This method finds the closest sample to the query point, and uses that closest sample as a starting point for the [Brent's Method](http://en.wikipedia.org/wiki/Brent%27s_method).

Spline Length Calculator
=============
The Spline Length calculator object is used to answer the question "How long is the spline segment between these two T values?"

To create a Spline Length Calculator, import `spline_library/splinelengthcalculator.h` and create a new `SplineLengthCalcuator` object by passing a `std::shared_ptr<Spline>` to the constructor. `SplineLengthCalculator lengthCalculator(myspline);`

### double findLength(double beginT, double endT, bool useShortestPath=false, double epsilon=0.0025) const
This method recursively subdivides the specified segment until the ssubdivisions stop improving the estimate of the length of the segment. It then returns the sum of the straight-line length of those subdivided segments.

If the spline is a non-looping spline, the useShortestPath parameter has no effect. If the spline is a looping spline and useShortestPath parameter is true, this method will compute the shortest path around the spline, either forwards or backwards, rather than blindly returning the distance from beginT to endT.

The epsilon sets the error tolerance. A value of 0 will make the calculation as precise as it can realistically be, but the computation will be slow value of 1 will provide a very rough, inacurate estimate, but will be very fast to compute. The default is meant to provide a reasonable mix of the two, being precise enough for most uses wile also being fast to compute.


