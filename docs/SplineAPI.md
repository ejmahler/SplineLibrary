Spline API
=============

The `Spline` class is a pure-virtual abstract class that serves as the base class of all spline types. This file documents its public methods. All interaction with splines is intended happen via this parent class.

#### Vector3D getPosition(t) const
This method computes the interpolated position at T.

The interpolated position is always continuous, IE there are never any "gaps" in the spline.

For looping splines, if the input T is outside of the allowable range (Ie less than zero or larger than the maximum T), the method will perform circular modular arithmetic to bring the T value into range. For example, if the input T value is -1, it will be adjusted to maxT - 1 before interpolation. This underlines the cyclic nature of looping splines.

For non-looping splines, if the input T is outside of the allowable range, the results are undefined.

#### Vector3D getTangent(t) const
This method interpolates the position and derivative of the position (aka the tangent) at T, and returns a struct containing both. There is no way to compute only the tangent.

The interpolated tangent is always continuous.

The derivative is determined analytically at compile-time, as opposed to numerically at run-time, so it is fast and numerically stable to compute.

The behavior when the T value is out of range is the same as for the getPosition method.

#### Vector3D getCurvature(t) const
This method interpolates the position, the first derivative of position (AKA the tangent), the second derivative of position (AKA the curvature), and returns a struct containing all three. There is no way to compute only the curvature.

The interpolated curvature is sometimes continuous - see the [Spline Types](SplineTypes.md) entry for your chosen spline type to find out whether or not its curvature is continuous.

The derivatives are determined analytically at compile-time, as opposed to numerically at run-time, so they are fast and numerically stable to compute.

The behavior when the T value is out of range is the same as for the getPosition method.

#### Vector3D getWiggle(t) const
This method interpolates the position, the first derivative of position (AKA the tangent), the second derivative of position (AKA the curvature), the third derivative of position (AKA the wiggle), and returns a struct containing all four. There is no way to compute only the wiggle.

For all current spline types, the wiggle is never continuous from segment to segment. For cubic splines, it is always a constant within each segment, although it may change from segment to segment.

The derivatives are determined analytically at compile-time, as opposed to numerically at run-time, so they are fast and numerically stable to compute.

The behavior when the T value is out of range is the same as for the getPosition method.

#### double getT(int index) const
The provided index should be an index into the original vector of control points provided to the constructor. The return value is the T value that the spline has assigned to that control point. This is very useful when using the [Centripetal Catmull-Rom Spline](SplineTypes.md#centripetal-catmull-rom-spline) or any other spline type where T values for points aren't necessarily evenly spaced.

The behavior when the index is out of range is undefined.

#### double getMaxT() const
This method returns the largest in-range T value.

All spline types normalize T values so that the maximum T value is equal to the index of the corresponding point.
* For non-looping splines, the maximum T depends on the spline type, but is generally equal to the number of points that were actually used, minus one. For example, a non-looping Catmull-Rom spline uses `size - 2` points internally - the other two are used only to calculate tangents. In this case, the maximum T value will be `size - 3`.
* For looping splines, the maximum T is always equal to the size of the input vector of points.

The minimum T value is always 0.

#### bool isLooping() const
Returns true if this spline is a looping spline, and false if this is a non-looping spline.
