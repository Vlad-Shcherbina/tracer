#ifndef __MATRIX__
#define __MATRIX__

#include <memory.h>
#include	"Vector3D.h"

class Matrix
{
public:
	float x [4][4];

	Matrix () {}
	Matrix ( float );
	Matrix ( const Matrix& m )
	{
		memcpy ( & x [0][0], &m.x [0][0], 16*sizeof ( float ) );
	}

	Matrix& operator += ( const Matrix& );
	Matrix& operator -= ( const Matrix& );
	Matrix& operator *= ( const Matrix& );
	Matrix& operator *= ( float );
	Matrix& operator /= ( float );

	void	invert ();
	void	transpose ();

	friend	Matrix   operator + ( const Matrix&, const Matrix& );
	friend	Matrix   operator - ( const Matrix&, const Matrix& );
	friend	Matrix   operator * ( const Matrix&, float );
	friend  Matrix   operator * ( float,         const Matrix& );
	friend	Matrix   operator * ( const Matrix&, const Matrix& );
	friend	Vector3D operator * ( const Matrix&, const Vector3D& );

	static Matrix	translate ( const Vector3D& );
	static Matrix	scale     ( const Vector3D& );
	static Matrix	rotateX   ( float );
	static Matrix	rotateY   ( float );
	static Matrix	rotateZ   ( float );
	static Matrix	rotate    ( const Vector3D& v, float );
	static Matrix	mirrorX   ();
	static Matrix	mirrorY   ();
	static Matrix	mirrorZ   ();
};


#endif
