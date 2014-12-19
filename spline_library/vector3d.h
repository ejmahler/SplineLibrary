#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>

class Vector3D
{
public:
	Vector3D();
	Vector3D(double x, double y, double z);

    inline double x() const;
    inline double y() const;
    inline double z() const;

    inline Vector3D &operator+=(const Vector3D &v);
    inline Vector3D &operator-=(const Vector3D &v);
    inline Vector3D &operator*=(double s);
    inline Vector3D &operator/=(double s);

	friend inline const Vector3D operator+(const Vector3D &left, const Vector3D &right);
    friend inline const Vector3D operator-(const Vector3D &left, const Vector3D &right);
    friend inline const Vector3D operator*(double s, const Vector3D &v);
    friend inline const Vector3D operator*(const Vector3D &v, double s);
    friend inline const Vector3D operator-(const Vector3D &v);
    friend inline const Vector3D operator/(const Vector3D &v, double s);

    friend inline bool operator==(const Vector3D &left, const Vector3D &right);
    friend inline bool operator!=(const Vector3D &left, const Vector3D &right);

    inline double length() const;
    inline double lengthSquared() const;

    inline Vector3D normalized() const;

    inline static double dotProduct(const Vector3D& left, const Vector3D& right);

private:
	double xp;
	double yp;
	double zp;
};

inline Vector3D::Vector3D()
	:xp(0.0), yp(0.0), zp(0.0)
{}

inline Vector3D::Vector3D(double xpos, double ypos, double zpos)
	:xp(xpos), yp(ypos), zp(zpos)
{}

inline double Vector3D::x() const { return xp; }
inline double Vector3D::y() const { return yp; }
inline double Vector3D::z() const { return zp; }

inline Vector3D &Vector3D::operator+=(const Vector3D &v)
{
    xp += v.xp;
    yp += v.yp;
    zp += v.zp;
    return *this;
}

inline Vector3D &Vector3D::operator-=(const Vector3D &v)
{
    xp -= v.xp;
    yp -= v.yp;
    zp -= v.zp;
    return *this;
}

inline Vector3D &Vector3D::operator*=(double s)
{
    xp *= s;
    yp *= s;
    zp *= s;
    return *this;
}

inline Vector3D &Vector3D::operator/=(double s)
{
    xp /= s;
    yp /= s;
    zp /= s;
    return *this;
}

inline const Vector3D operator+(const Vector3D &left, const Vector3D &right)
{
    return Vector3D(left.xp + right.xp, left.yp + right.yp, left.zp + right.zp);
}

inline const Vector3D operator-(const Vector3D &left, const Vector3D &right)
{
    return Vector3D(left.xp - right.xp, left.yp - right.yp, left.zp - right.zp);
}

inline const Vector3D operator*(double s, const Vector3D &v)
{
    return Vector3D(v.xp * s, v.yp * s, v.zp * s);
}

inline const Vector3D operator*(const Vector3D &v, double s)
{
    return Vector3D(v.xp * s, v.yp * s, v.zp * s);
}

inline const Vector3D operator-(const Vector3D &v)
{
    return Vector3D(-v.xp, -v.yp, -v.zp);
}

inline const Vector3D operator/(const Vector3D &v, double s)
{
    return Vector3D(v.xp / s, v.yp / s, v.zp / s);
}

inline bool operator==(const Vector3D &left, const Vector3D &right)
{
	return left.xp == right.xp
		&& left.yp == right.yp
		&& left.zp == right.zp;
}

inline bool operator!=(const Vector3D &left, const Vector3D &right)
{
	return left.xp != right.xp
		|| left.yp != right.yp
		|| left.zp != right.zp;
}

inline Vector3D Vector3D::normalized() const
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

inline double Vector3D::dotProduct(const Vector3D& v1, const Vector3D& v2)
{
    return v1.xp * v2.xp + v1.yp * v2.yp + v1.zp * v2.zp;
}

inline double Vector3D::length() const
{
    return sqrt(xp * xp + yp * yp + zp * zp);
}

inline double Vector3D::lengthSquared() const
{
    return xp * xp + yp * yp + zp * zp;
}


#endif // VECTOR3D_H
