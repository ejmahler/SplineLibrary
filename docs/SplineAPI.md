Spline API
=============

The `Spline` class is an abstract class that serves as the base class of all spline types. This file documents its public methods.

The Spline class (and all subclasses) is a template class. It takes two template arguments: The first is the type of data to be interpolated (IE, a 2d vector, or a 3d vector). It expects the vector to have an API matching Qt's QVector2D class, but will support data with any dimension. Custom classes may be used as long as they provide the same public methods as QVector2D. Custom classes can ven use double instead of float as the internal data type if more precision is needed.

The second templae marameter is the floating point type to use for internal calculations (IE, float, double, some BigDecimal class). This defaults to float and it's safe to leave this a float for most applications, but if you use an InterpolationType custom class that stores its data as doubles, you'll get more precision by telling the Spline to use doubles as well.

#### getPosition(t)
This method computes the interpolated position at T.

The interpolated position is always continuous, IE there are never any "gaps" in the spline.

For looping splines, if the input T is outside of the allowable range (Ie less than zero or larger than the maximum T), the method will perform circular modular arithmetic to bring the T value into range. For example, if the input T value is -1, it will be adjusted to maxT - 1 before interpolation. This underlines the cyclic nature of looping splines.

For non-looping splines, if the input T is outside of the allowable range, the results are undefined.

#### getTangent(t)
This method interpolates the position and derivative of the position (aka the tangent) at T, and returns a struct containing both. While it's possible to compute only the tangent, the API doesn't support it, in the interest of simplifying the API.

The interpolated tangent is always continuous.

The derivative is determined analytically at compile-time, as opposed to numerically at run-time, so it is fast and numerically stable to compute.

The behavior when the T value is out of range is the same as for the getPosition method.

#### getCurvature(t)
This method interpolates the position, the first derivative of position (AKA the tangent), the second derivative of position (AKA the curvature), and returns a struct containing all three. While it's possible to compute only the curvature, the API doesn't support it, in the interest of simplifying the API.

The interpolated curvature is sometimes continuous - see the [Spline Types](SplineTypes.md) entry for your chosen spline type to find out whether or not its curvature is continuous.

The derivatives are determined analytically at compile-time, as opposed to numerically at run-time, so they are fast and numerically stable to compute.

The behavior when the T value is out of range is the same as for the getPosition method.

#### getWiggle(t) const
This method interpolates the position, the first derivative of position (AKA the tangent), the second derivative of position (AKA the curvature), the third derivative of position (AKA the wiggle), and returns a struct containing all four. While it's possible to compute only the wiggle, the API doesn't support it, in the interest of simplifying the API.

For all current spline types, the wiggle is never continuous from segment to segment. For cubic splines, it is always a constant within each segment, although it may change from segment to segment.

The derivatives are determined analytically at compile-time, as opposed to numerically at run-time, so they are fast and numerically stable to compute.

The behavior when the T value is out of range is the same as for the getPosition method.

#### arcLength(a, b) const
This method computes the arc length between a and b. IE, if you traceda path with your finger along the spline from a to b, how much distance would it cover?

This is found by numerically computing the integral of the magnitude of the tangent. In real world terms, it computes the tangent at several points between a and b and then combines the results.

For looping splines, it will use modular arithmetic to ensure that a and b are less than one "circuit" away from each other. Notably, this means that `arcLength(0, maxT)` will return 0 for looping splines, because it detects that 0 to maxT is a complete circuit and removes it. If you want to compute the length of the whole spline, use `totalLength()` instead.

#### totalLength() const
This method computes the arc length of the entire spline, from beginning to end. IE, if you traceda path with your finger along the spline from start to end, how much distance would it cover?

This is found by numerically computing the integral of the magnitude of the tangent. In real world terms, it computes the tangent at several points between the start and end and then combines the results.

For computing the total length of non-looping splines, calling `totalLength()` is preferred over calling `arcLength(0, maxT)` because it's slightly faster.

#### getT(int index) const
The provided index should be an index into the original vector of control points provided to the constructor. The return value is the T value that the spline has assigned to that control point. This is very useful when using the [Centripetal Catmull-Rom Spline](SplineTypes.md#centripetal-catmull-rom-spline) or any other spline type where T values for points aren't necessarily evenly spaced.

The behavior when the index is out of range is undefined.

#### getMaxT() const
This method returns the largest in-range T value.

All spline types normalize T values so that the maximum T value is equal to the index of the corresponding point.
* For non-looping splines, the maximum T depends on the spline type, but is generally equal to the number of points that were actually used, minus one. For example, a non-looping Catmull-Rom spline uses `size - 2` points internally - the other two are used only to calculate tangents. In this case, the maximum T value will be `size - 3`.
* For looping splines, the maximum T is always equal to the size of the input vector of points.

The minimum T value is always 0.

#### isLooping() const
Returns true if this spline is a looping spline, and false if this is a non-looping spline.
