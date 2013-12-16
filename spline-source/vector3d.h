#ifndef VECTOR3D_H
#define VECTOR3D_H

class Vector3D
{
public:
	Vector3D();
	Vector3D(double x, double y, double z);

	double x() const;
	double y() const;
	double z() const;

	Vector3D &operator+=(const Vector3D &v);
    Vector3D &operator-=(const Vector3D &v);
    Vector3D &operator*=(double s);
    Vector3D &operator/=(double s);

	friend inline const Vector3D operator+(const Vector3D &left, const Vector3D &right);
    friend inline const Vector3D operator-(const Vector3D &left, const Vector3D &right);
    friend inline const Vector3D operator*(double s, const Vector3D &v);
    friend inline const Vector3D operator*(const Vector3D &v, double s);
    friend inline const Vector3D operator-(const Vector3D &v);
    friend inline const Vector3D operator/(const Vector3D &v, double s);

    friend inline bool operator==(const Vector3D &left, const Vector3D &right);
    friend inline bool operator!=(const Vector3D &left, const Vector3D &right);

	double length() const;
    double lengthSquared() const;

    Vector3D normalized() const;

    static double dotProduct(const Vector3D& left, const Vector3D& right);

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


#endif // VECTOR3D_H
