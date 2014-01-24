#include <math.h>
#include "Matrix.h"

Matrix :: Matrix ( float v )
{
	for ( int i = 0; i < 4; i++)
		for ( int j = 0; j < 4; j++)
			x [i][j] = (i == j) ? v : 0.0;

	x [3][3] = 1;
}

void	Matrix :: invert ()
{
	Matrix	out ( 1 );

	for ( int i = 0; i < 4; i++ )
	{
		float	d = x [i][i];

		if ( d != 1.0)
		{
			for ( int j = 0; j < 4; j++ )
			{
				out.x [i][j] /= d;
				x [i][j]     /= d;
			}
		}

		for ( int j = 0; j < 4; j++ )
		{
			if ( j != i )
			{
				if ( x [j][i] != 0.0)
				{
					float	mulBy = x[j][i];

					for ( int k = 0; k < 4; k++ )
					{
						x [j][k]     -= mulBy * x [i][k];
						out.x [j][k] -= mulBy * out.x [i][k];
					}
				}
			}
		}
	}

	*this = out;
}

void	Matrix :: transpose ()
{
	float	t;

	for ( register int i = 0; i < 4; i++ )
		for ( register int j = i; j < 4; j++ )
			if ( i != j )
			{
				t        = x [i][j];
				x [i][j] = x [j][i];
				x [j][i] = t;
			}
}

Matrix&	Matrix :: operator += ( const Matrix& a )
{
	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
			x [i][j] += a.x [i][j];

	return *this;
}

Matrix&	Matrix :: operator -= ( const Matrix& a )
{
	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
			x [i][j] -= a.x [i][j];

	return *this;
}

Matrix&	Matrix :: operator *= ( float v )
{
	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
			x [i][j] *= v;

	return *this;
}

Matrix&	Matrix :: operator *= ( const Matrix& a )
{
	Matrix	res ( *this );

	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
		{
			float sum = 0;

			for ( int k = 0; k < 4; k++ )
				sum += res.x [i][k] * a.x [k][j];

			x [i][j] = sum;
		}

    return *this;
}

Matrix	operator + ( const Matrix& a, const Matrix& b )
{
	Matrix	res;

	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
			res.x [i][j] = a.x [i][j] + b.x [i][j];

	return res;
}

Matrix	operator - ( const Matrix& a, const Matrix& b )
{
	Matrix	res;

	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
			res.x [i][j] = a.x [i][j] - b.x [i][j];

	return res;
}

Matrix	operator * ( const Matrix& a, const Matrix& b )
{
	Matrix	res;

	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
		{
			float sum = 0;

			for ( register int k = 0; k < 4; k++ )
				sum += a.x [i][k] * b.x [k][j];

			res.x [i][j] = sum;
		}

	return res;
}

Matrix	operator * ( const Matrix& a, float v )
{
	Matrix	res;

	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
			res.x [i][j] = a.x [i][j] * v;

	return res;
}

Matrix	operator * ( float v, const Matrix& a )
{
	Matrix	res;

	for ( register int i = 0; i < 4; i++ )
		for ( register int j = 0; j < 4; j++ )
			res.x [i][j] = a.x [i][j] * v;

	return res;
}

Vector3D operator * ( const Matrix& m, const Vector3D& v )
{
	Vector3D res;

	res.x = m.x [0][0] * v.x + m.x [0][1] * v.y + m.x [0][2] * v.z + m.x [0][3];
	res.y = m.x [1][0] * v.x + m.x [1][1] * v.y + m.x [1][2] * v.z + m.x [1][3];
	res.z = m.x [2][0] * v.x + m.x [2][1] * v.y + m.x [2][2] * v.z + m.x [2][3];

	float	denom = m.x [3][0] * v.x + m.x [3][1] * v.y +  m.x [3][2] * v.z + m.x [3][3];

	if ( denom != 1.0 )
		res /= denom;

	return res;
}

//////////////////////// Derived functions /////////////////////////////

Matrix	Matrix::translate ( const Vector3D& loc )
{
	Matrix	res ( 1 );

	res.x [0][3] = loc.x;
	res.x [1][3] = loc.y;
	res.x [2][3] = loc.z;

	return res;
}

Matrix	Matrix::scale ( const Vector3D& v )
{
	Matrix	res ( 1 );

	res.x [0][0] = v.x;
	res.x [1][1] = v.y;
	res.x [2][2] = v.z;

	return res;
}

Matrix	Matrix::rotateX ( float angle )
{
	Matrix res ( 1 );
	float  cosine = cos ( angle );
	float  sine   = sin ( angle );

	res.x [1][1] = cosine;
	res.x [1][2] = -sine;
	res.x [2][1] = sine;
	res.x [2][2] = cosine;

	return res;
}

Matrix	Matrix::rotateY ( float angle )
{
	Matrix res ( 1 );
	float  cosine = cos ( angle );
	float  sine   = sin ( angle );

	res.x [0][0] = cosine;
	res.x [0][2] = -sine;
	res.x [2][0] = sine;
	res.x [2][2] = cosine;

	return res;
}

Matrix	Matrix::rotateZ ( float angle )
{
	Matrix res ( 1 );
	float  cosine = cos ( angle );
	float  sine   = sin ( angle );

	res.x [0][0] = cosine;
	res.x [0][1] = -sine;
	res.x [1][0] = sine;
	res.x [1][1] = cosine;

	return res;
}

Matrix	Matrix::rotate ( const Vector3D& axis1, float angle )
{
	Matrix res ( 1 );
	Vector3D axis=axis1; axis.normalize();
	
	float  cosine = cos ( angle );
	float  sine   = sin ( angle );

	res.x [0][0] = axis.x * axis.x + ( 1 - axis.x * axis.x ) * cosine;
	res.x [1][0] = axis.x * axis.y * ( 1 - cosine ) + axis.z * sine;
	res.x [2][0] = axis.x * axis.z * ( 1 - cosine ) - axis.y * sine;
	res.x [3][0] = 0;

	res.x [0][1] = axis.x * axis.y * ( 1 - cosine ) - axis.z * sine;
	res.x [1][1] = axis.y * axis.y + ( 1 - axis.y * axis.y ) * cosine;
	res.x [2][1] = axis.y * axis.z * ( 1 - cosine ) + axis.x * sine;
	res.x [3][1] = 0;

	res.x [0][2] = axis.x * axis.z * ( 1 - cosine ) + axis.y * sine;
	res.x [1][2] = axis.y * axis.z * ( 1 - cosine ) - axis.x * sine;
	res.x [2][2] = axis.z * axis.z + ( 1 - axis.z * axis.z ) * cosine;
	res.x [3][2] = 0;

	res.x [0][3] = 0;
	res.x [1][3] = 0;
	res.x [2][3] = 0;
	res.x [3][3] = 1;

	return res;
}

Matrix	Matrix::mirrorX ()
{
	Matrix	res ( 1 );

	res.x [0][0] = -1;

	return res;
}

Matrix	Matrix::mirrorY ()
{
	Matrix	res ( 1 );

	res.x [1][1] = -1;

	return res;
}

Matrix	Matrix::mirrorZ ()
{
	Matrix	res ( 1 );

	res.x [2][2] = -1;

	return res;
}
