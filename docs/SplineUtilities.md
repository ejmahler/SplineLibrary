This file outlines some utility classes for interacting with splines.

Spline Inverter
=============
The Spline Inverter object is used to answer the question "Given a data point (called the Query Point), what t value brings the spline closest to that point?" This can be used for query points which are on the spline itself, but it can also be used for points that are near the spline, or even for query points that are far away.

An example use case is an AI system for a race track. The spline could present the optimal line through a curve in the track, and the AI system could use its current position as a query point to the Spline Inverter to get the closest point on the spline. This will tell the AI where it should be via the interpolated position, what direction it should be going via the interpolated tangent tangent, if it should be turning via the interpolated curvature, etc.

This must be computed numerically - there is no reasonable analytic solution like there is for derivatives.

To create a Spline Inverter, import `spline_library/utils/splineinverter.h`, and create a new `SplineInverter` object by passing a reference to a Spline to the constructor.
```c++
std::vector<QVector2D> splinePoints = ...;
UniformCRSpline<QVector2D> mySpline(splinePoints);
SplineInverter<QVector2D> inverter(mySpline);
```

Alternatively, if you have a shared_ptr to a Spline, you can do the following:
```c++
std::shared_ptr<Spline<QVector2D>> mySplinePtr = ...;
SplineInverter<QVector2D> inverter(*mySplinePtr);
```

The SplineInverter stores a reference to the spline, so it should not live longer than the spline it refers to.

In the SplineInverter constructor, it takes "samples" of the spline at regular intervals. By default it takes 10 samples per T, but this can be changed via a constructor parameter. When given a query point, it first finds the closest sample to the query point, then uses that sample location as the starting point for a refining algorithm.

### findClosestT(queryPoint) const
This method finds the closest sample to the query point, and uses that closest sample as a starting point for [Brent's Method](http://en.wikipedia.org/wiki/Brent%27s_method).

Example:
```c++
SplineInverter<QVector2D> inverter = ...;
float t = inverter.findClosestT(QVector2D(5, 1));
```


Arc Length Solver
=============
The arc length solver methods, found in `spline_library/utils/arclength.h` all deal with a similar question: Given a starting t value on the spline and a desired arc length, what secondary T value will yield my desired arc length? All methods listed here will accept any spline type. They will accept references to the parent Spline class, but they're all template functions on spline type, so it's possible to avoid virtual function calls by passing in a reference to a concrete spline type.

### ArcLength::solveLength(const spline&, a, desiredLength)
Given a spline, a starting T value on the spline `a`, and a desired arc length `desiredLength`, compute and return `b` that satisfies the condition `spline.arcLength(a, b) ~= desiredLength`

If the desired length goes past the end of the spline (IE if b would be greater than maxT), maxT is returned.

Example:
```c++
std::vector<QVector2D> splinePoints = ...;
UniformCRSpline<QVector2D> mySpline(splinePoints);

float a = 1.2f;
float desiredArcLength = 1.5f;

float b = ArcLength::solveLength(mySpline, a, desiredArcLength);
```


### ArcLength::partition(const spline&, desiredLength)
Given a spline and a desired arc length `desiredLength`, partition the spline into as many pieces as possible, such that each piece has arc length `desiredLength`. Returns a `std::vector` of T values marking the beginning/end of each partitioned piece.

IE, if it returns `std::vector<float>{0.0, 1.2, 2.6}` this signifies that the spline was partitioned into two pieces, with the first starting at T = 0 and ending at T = 1.2, and the second beginning at T = 1.2 and ending at T = 2.6.

If `desiredLength` does not cleanly divide the spline's total arc length (it almost certainly won't due to floating point inaccuracy), the "remainder" is not included in the list, so the final value will be the T value at the end of the final full-length piece. If you want a guarantee that the entire spline is partitioned and there is no remainder, see `ArcLength::partitionN` below.

Example:
```c++
std::vector<QVector2D> splinePoints = ...;
UniformCRSpline<QVector2D> mySpline(splinePoints);

float desiredArcLength = 1.5f;

std::vector<float> partitionBoundaries = ArcLength::partition(mySpline, desiredArcLength);
```


### ArcLength::partitionN(const spline&, n)
Given a spline and a number pf pieces `n`, divide the spline into `n` pieces such that each piece has the same arc length. Returns a `std::vector` containing n+1 T values marking the beginning/end of each partitoned piece.

The first element in the returned vector is always 0, and the last element is always maxT.

Example:
```c++
std::vector<QVector2D> splinePoints = ...;
UniformCRSpline<QVector2D> mySpline(splinePoints);

size_t n = 6;

std::vector<float> partitionBoundaries = ArcLength::partitionN(mySpline, n);
```
