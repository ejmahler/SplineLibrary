This file outlines some utility classes for interacting with splines.

Spline Inverter
=============
The Spline Inverter object is used to answer the question "Given a data point (called the Query Point), what t value brings the spline closest to that point?" This can be used for query points which are on the spline itself, but it can also be used for points that are near the spline, or even for query points that are far away.

An example use case is an AI system for a race track. The spline could present the optimal line through a curve in the track, and the AI system could use its current position as a query point to the Spline Inverter to get the closest point on the spline. This will tell the AI where it should be, what direction it should be going via the tangent, how the optimal line is turning via the curvature, etc.

This must be computed numerically - there is no reasonable analytic solution like there is for derivatives.

To create a Spline Inverter, import `spline_library/splineinverter.h`, and create a new `SplineInverter` object by passing a reference to a Spline to the constructor. Example: `SplineInverter inverter(*mySpline.get());` where mySpline is a `std::shared_ptr<Spline>`.

In the SplineInverter constructor, it takes "samples" of the spline at regular intervals. By default it takes 10 samples per T, but this can be changed via a constructor parameter. When given a query point, it first finds the closest sample to the query point, then uses that sample location as the starting point for a refining algorithm.

### findClosestT(queryPoint) const
This method finds the closest sample to the query point, and uses that closest sample as a starting point for [Brent's Method](http://en.wikipedia.org/wiki/Brent%27s_method).