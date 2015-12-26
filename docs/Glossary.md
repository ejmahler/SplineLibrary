Glossary
=============

This page contains definitions for terms that are commonly used in the documentation.

### Spline
A spline is a piecewise function used to smoothly interpolate through or approximate a data set. A spline takes an input value named `t` and returns an interpolated value. The input t values typically correlate with the original data set, IE an input t value equal to n will typically return the n-th input point, and an input t value of n + 0.5 will return a value somewhere in the vicinity of the mid point between the n-th point and the (n+1)-th point.

A core goal of a spline is not just interpolation, but smooth interpolation. This means that, rather than the spline sharply changing direction, it will curve around from one point to the next.

### Piecewise
A spline is said to be piecewise because it is not composed of one single mathematical formula. Instead, we say "if the input t is between 0 and 1, use one formula. if the input t is between 1 and 2, use a different formula. if the input t is between 2 and 3, use a still different formula" and so on for every set of input points. The differences in spline types come from how they construct each of these formulas.

### Derivative
The derivative comes from calculus - it describes the rate of change of a function. It answers the question "Given a mathematical formula, how does that formula's output change as the input changes?"

In the context of a spline, the first derivative answers the question "in what direction is the spline moving at `t`?". We call this the "tangent".

We can also take the derivative of a derivative. This a "second derivative" of the original function. It describes the rate of change of the first derivative. In the context of a spline, this answers the questions "Is the spline speeding up or slowing down at `t`? Is it turning left or turning right?" and is called the "curvature".

### Continuous
In layman's terms, a function is said to be "continuous" if there are no "breaks" or "jumps" in the function. It is always possible to trace a path with our finger from the start of a spline to the end of a spline without lifting our finger.

Derivatives of a function can be continuous or non-continous as well. We can visually check if a spline's derivative is continuous in the demo. If a spline's derivative wasn't continuous, the spline would appear to sharply change direction, rather than forming a smooth curve.

As said before, a core goal of a spline is not just interpolation, but smooth interpolation. In order to meat this goal, a spline must be continuous, and it must have a continuous first derivative.

### Continuous Curvature
Second derivatives can also be continuous or non-continous, but this is not a priority for many spline types. It is much harder to visually verify that a spline's curvature is continuous, but it can be a useful property for some applications, so the [Spline Types](SplineTypes.md) page lists whether or not the curvature is continuous for each type.

### Local Control
Local control is an optional property of splines that is very desirable for many applications. For a spline with local control, moving a data point will only affect the spline segments close to that point, leaving all others unaffected. This can be seen in the demo when looking at a `CRSpline`, which has local control: If we move one data point, two segments to the left and two segments to the right will move, but all other parts of the spline will remain untouched.

Conversely, the `NaturalSpline` is an example of a spline that does not have local control. If we move a data point in the demo, every single spline segment is affected, rather than just the ones near the data point.

It is interesting to point out here that the Natural Spline has a continuous curvature, while the Catmull-Rom Spline does not, this highlighting the fact that choosing a spline type is about balancing tradeoffs between different features of each type. 
