//
// Basic class for 3d vectors.
// Supports al basic vector operations
//
// Author: Alex V. Boreskoff
//
// Last change: 09/12/1999
//

#ifndef	__VECTOR3D__
#define	__VECTOR3D__

//#define FAST_VECTORS

#include	<math.h>

class	Vector3D
{
public:
	float	x, y, z;

	Vector3D () {}
	Vector3D ( float px, float py, float pz )
	{
		x = px;
		y = py;
		z = pz;
	}

#ifdef FAST_VECTORS
private:
#endif
	Vector3D ( const Vector3D& v )
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}
public:

	Vector3D& operator = ( const Vector3D& v )
	{
		x = v.x;
		y = v.y;
		z = v.z;

		return *this;
	}

	Vector3D operator + () const
	{
		return *this;
	}

	Vector3D operator - () const
	{
		return Vector3D ( -x, -y, -z );
	}

	Vector3D& operator += ( const Vector3D& v )
	{
		x += v.x;
		y += v.y;
		z += v.z;

		return *this;
	}

	Vector3D& operator -= ( const Vector3D& v )
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;

		return *this;
	}

	Vector3D& operator *= ( const Vector3D& v )
	{
		x *= v.x;
		y *= v.y;
		z *= v.z;

		return *this;
	}

	Vector3D& operator *= ( float f )
	{
		x *= f;
		y *= f;
		z *= f;

		return *this;
	}

	Vector3D& operator /= ( const Vector3D& v )
	{
		x /= v.x;
		y /= v.y;
		z /= v.z;

		return *this;
	}

	Vector3D& operator /= ( float f )
	{
		x /= f;
		y /= f;
		z /= f;

		return *this;
	}

	float& operator [] ( int index )
	{
		return * ( index + &x );
	}

	int	operator == ( const Vector3D& v ) const
	{
		return x == v.x && y == v.y && z == v.z;
	}

	int	operator != ( const Vector3D& v ) const
	{
		return x != v.x || y != v.y || z != v.z;
	}

	operator float * ()
	{
		return &x;
	}

	float	length () const
	{
		return (float) sqrt ( x * x + y * y + z * z );
	}

	Vector3D&	normalize ()
	{
		return (*this) /= length ();
	}

    float   maxLength () 
    {
        return max3 ( fabs (x), fabs (y), fabs (z) );
    }

	friend Vector3D operator + ( const Vector3D&, const Vector3D& );
	friend Vector3D operator - ( const Vector3D&, const Vector3D& );
	friend Vector3D operator * ( const Vector3D&, const Vector3D& );
	friend Vector3D operator * ( float,           const Vector3D& );
	friend Vector3D operator * ( const Vector3D&, float );
	friend Vector3D operator / ( const Vector3D&, float );
	friend Vector3D operator / ( const Vector3D&, const Vector3D& );
	friend float    operator & ( const Vector3D&, const Vector3D& );
	friend Vector3D operator ^ ( const Vector3D&, const Vector3D& );

private:
    float   max3 ( float a, float b, float c )
    {
        return a > b ? (a > c ? a : (b > c ? b : c)) :
                       (b > c ? b : (a > c ? a : c));
    }
};

inline Vector3D operator + ( const Vector3D& u, const Vector3D& v )
{
	return Vector3D ( u.x + v.x, u.y + v.y, u.z + v.z );
}

inline Vector3D operator - ( const Vector3D& u, const Vector3D& v )
{
	return Vector3D ( u.x - v.x, u.y - v.y, u.z - v.z );
}

inline Vector3D operator * ( const Vector3D& u, const Vector3D& v )
{
	return Vector3D ( u.x*v.x, u.y*v.y, u.z * v.z );
}

inline Vector3D operator * ( const Vector3D& v, float a )
{
	return Vector3D ( v.x*a, v.y*a, v.z*a );
}

inline Vector3D operator * ( float a, const Vector3D& v )
{
	return Vector3D ( v.x*a, v.y*a, v.z*a );
}

inline Vector3D operator / ( const Vector3D& u, const Vector3D& v )
{
	return Vector3D ( u.x/v.x, u.y/v.y, u.z/v.z );
}

inline Vector3D operator / ( const Vector3D& v, float a )
{
	return Vector3D ( v.x/a, v.y/a, v.z/a );
}

inline float operator & ( const Vector3D& u, const Vector3D& v )
{
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline Vector3D operator ^ ( const Vector3D& u, const Vector3D& v )
{
	return Vector3D (u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}

#endif
