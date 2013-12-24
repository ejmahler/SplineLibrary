#include "vector3d.h"

#include <cmath>

Vector3D Vector3D::normalized() const
{
    // Need some extra precision if the length is very small.
    double lengthSquared =	xp * xp +
							yp * yp +
							zp * zp;
    if (lengthSquared < .00000001)
		return Vector3D();
	else
	{
		double invLength = 1 / sqrt(lengthSquared);
        return Vector3D(xp * invLength, yp * invLength, zp * invLength);
	}
}

double Vector3D::dotProduct(const Vector3D& v1, const Vector3D& v2)
{
    return v1.xp * v2.xp + v1.yp * v2.yp + v1.zp * v2.zp;
}

double Vector3D::length() const
{
    return sqrt(xp * xp + yp * yp + zp * zp);
}

double Vector3D::lengthSquared() const
{
    return xp * xp + yp * yp + zp * zp;
}