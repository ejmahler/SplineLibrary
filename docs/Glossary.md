Glossary
=============

This page contains definitions for terms that are commonly used in the documentation.

### Spline
A spline is a piecewise function used to smoothly interpolate through or approximate a data set. A spline takes an input value named `t` and returns an interpolated value. The input t values typically correlate with the original data set, IE an input t value equal to n will typically return the n-th input point, and an input t value of n + 0.5 will return a value somewhere in the vicinity of the mid point between the n-th point and the (n+1)-th point.

A core goal of a spline is not just interpolation, but smooth interpolation. This means that, rather than the spline sharply changing direction, it will curve around from one point to the next.

### Piecewise
A spline is said to be piecewise because it is not composed of one single mathematical formula. Instead, we say "if the input t is between 0 and 1, use this formula. if the input t is between 1 and 2, use this formula" and so on for every set of input points. The differences in spline types come from how they construct each of these formulas.

### Derivative
The derivative comes from calculus - it describes the rate of change of a function. It answers the question "Given a mathematical formula, how does that formula's output change?"

In the context of a spline, the first derivative answers the question "in what direction is the spline moving at `t`?". We call this the "tangent".

The derivative of a function is itself a function, so we can logically extend the idea of a derivative to this function too. When we take the derivative of a derivative, we obtain a "second derivative" - that is, the rate of change of the first derivative. In the context of a spline, this answers the questions "Is the spline speeding up or slowing down at `t`? Is it turning left or turning right?" and is called the "curvature".

### Continuous
In layman's terms, a function is said to be "continuous" if there are no "breaks" or "jumps" in the function. An informal mathematical definition might be `f is continuous if and only if, for all x, the limit of f(x) - f(x + h) as h goes to 0 is equal to 0`.

In the context of the spline, we can easily verify that the spline is continous, because it doesn't stop at one point and suddenly start up again somewhere else, with a gap in between. It is always possible to trace a path with our finger from the start of a spline to the end of a spline without lifting our finger.

Derivatives of a function can be continuous or non-continous as well. As said before, a core goal of a spline is not just interpolation, but smooth interpolation. Now that we've defined the derivative, we can more rigorously define smooth "smooth interpolation" - a function is smooth if both it and its derivative are continuous. Thus, all splines must be continuous, and all splines must have a continuous first derivative, if they are to meet the goal of smooth interpolation.

We can verify that a spline's derivative is continuous in the demo. If a spline's derivative wasn't continuous, the spline would appear to sharply change direction, rather than forming a smooth curve.

### Continuous Curvature
Second derivatives can also be continuous or non-continous, but this is not a priority for many spline types. It is much harder to visually verify that a spline's curvature is continuous, but it can be a useful property for some applications, so the [Spline Types](SplineTypes.md) page lists whether or not the curvature is continuous for each type.

### Local Control
Local control is an optional property of splines that is very desirable for many applications. For a spline with local control, moving a data point will only affect the spline segments near to that point, leaving all others affected. This can be seen in the demo when looking at a `CRSpline`, which has local control: If we move one data point, two segments to the left and to segments to the right will move, but all other parts of the spline will remain untouched.

Conversely, the `NaturalSpline` is an example of a spline that does not have local control. If we move a data point in the demo, every single spline segment is affected, rather than just the ones near the data point.

It is useful to point out here that the Natural Spline has a continuous curvature, while the Catmull-Rom Spline does not, this highlighting the fact that choosing a spline type is about balancing tradeoffs between different features of each type. 
