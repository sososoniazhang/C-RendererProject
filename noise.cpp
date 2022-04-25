#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <math.h>
#include "noise.h"
#include <algorithm>
#define PERLIN 1
#define VALUE 2
#define SIMPLEX 3
#define NORMAL_NOISE 15
#define SUM_NOISE 16
#define SUM_NOISE_ABS 17
#define SIN_NOISE 18
//#define MOD3 vec3(.1031,.11369,.13787)
using namespace std;

Noise::Noise() {};

Noise::~Noise() {};

/*Function for dot product*/
float dotProduct(Vec3f v1, Vec3f v2)
{
	float product = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

	return product;
}

/*Function for dot product*/
float dotProduct4(Vec4f v1, Vec4f v2)
{
	float product = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;

	return product;
}

float fract(float number)
{
	return number - floor(number);
}

Vec3f step(Vec3f edge, Vec3f x_point)
{
	Vec3f result;
	if (x_point.x < edge.x) {
		result.x = 0.0;
	}
	else {
		result.x = 1.0;
	}

	if (x_point.y < edge.y) {
		result.y = 0.0;
	}
	else {
		result.y = 1.0;
	}
	if (x_point.z < edge.z) {
		result.z = 0.0;
	}
	else {
		result.z = 1.0;
	}

	return result;
}

Vec3f vector_time_matrix(Vec3f vec, float mat[3][3], float cof)
{
	Vec3f row1 = { cof * mat[0][0], cof * mat[0][1], cof * mat[0][2] };
	Vec3f row2 = { cof * mat[1][0], cof * mat[1][1], cof * mat[1][2] };
	Vec3f row3 = { cof * mat[2][0], cof * mat[2][1], cof * mat[2][2] };

	Vec3f u;
	u.x = dotProduct(vec, row1);
	u.y = dotProduct(vec, row2);
	u.z = dotProduct(vec, row3);

	return u;
}


Vec3f Noise::random3(Vec3f c)
{
	Vec3f p1 = { 17.0, 59.40, 15.0 };
	float j = sinf(dotProduct(c, p1)/100.0);
	Vec3f result;
	result.z = (fract(512.0 * j) - 0.5);
	j *= 12.5;
	result.x = (fract(512.0 * j) - 0.5);
	j *= 12.5;
	result.y = (fract(512.0 * j) - 0.5);
	return result;
}

float Noise::simplex3d(Vec3f p)
{	 /* 1. find current tetrahedron T and it's four vertices */
	 /* s, s+i1, s+i2, s+1.0 - absolute skewed (integer) coordinates of T vertices */
	 /* x, x1, x2, x3 - unskewed coordinates of p relative to each of T vertices*/

	 /* calculate s and x */
	Vec3f F1 = { F3, F3, F3 };
	float num = dotProduct(p, F1);
	Vec3f s = { floor(p.x + num), floor(p.y + num), floor(p.z + num) };

	Vec3f G1 = { G3, G3, G3 };
	float num2 = dotProduct(s, G1);
	Vec3f x1 = { p.x - s.x + num2, p.y - s.y + num2, p.z - s.z + num2 };
	/* calculate i1 and i2 */
	Vec3f yzx = { x1.y, x1.z, x1.x };
	Vec3f x_yzx = { x1.x - yzx.x, x1.y - yzx.y, x1.z - yzx.z };
	Vec3f zero = { 0.0, 0.0, 0.0 };
	Vec3f e = step(zero, x_yzx);

	Vec3f zxy = { e.z, e.x, e.y };// vector e.z, e.x, e.y
	Vec3f one_zxy = { 1.0 - zxy.x, 1.0 - zxy.y, 1.0 - zxy.z };
	Vec3f i1 = { e.x * one_zxy.x, e.y * one_zxy.y, e.z * one_zxy.z };

	Vec3f one_e = { 1.0 - e.x, 1.0 - e.y, 1.0 - e.z };
	Vec3f zxy_oneE = { zxy.x * one_e.x, zxy.y * one_e.y, zxy.z * one_e.z };
	Vec3f i2 = { 1.0 - zxy_oneE.x,1.0 - zxy_oneE.y,1.0 - zxy_oneE.z };

	Vec3f x11 = { x1.x - i1.x + G3, x1.y - i1.y + G3 ,x1.z - i1.z + G3 };
	Vec3f x22 = { x1.x - i2.x + 2.0 * G3, x1.y - i2.y + 2.0 * G3 ,x1.z - i2.z + 2.0 * G3 };
	Vec3f x33 = { x1.x - 1.0 + 3.0 * G3,x1.y - 1.0 + 3.0 * G3,x1.z - 1.0 + 3.0 * G3 };
	/* 2. find four surflets and store them in d */
	Vec4f w;
	Vec4f d;
	/* calculate surflet weights */
	 /* w fades from 0.6 at the center of the surflet to 0.0 at the margin */
	w.x = max(0.6 - dotProduct(x1, x1), 0.0);
	w.y = max(0.6 - dotProduct(x11, x11), 0.0);
	w.z = max(0.6 - dotProduct(x22, x22), 0.0);
	w.w = max(0.6 - dotProduct(x33, x33), 0.0);
	/* calculate surflet components */
	/* calculate surflet components */
	Vec3f s_i1 = { s.x + i1.x, s.y + i1.y, s.z + i1.z };
	Vec3f s_i2 = { s.x + i2.x, s.y + i2.y, s.z + i2.z };
	Vec3f s_one = { s.x + 1.0, s.y + 1.0, s.z + 1.0 };
	d.x = dotProduct(random3(s), x1);
	d.y = dotProduct(random3(s_i1), x11);
	d.z = dotProduct(random3(s_i2), x22);
	d.w = dotProduct(random3(s_one), x33);

	/* multiply d by w^4 */
	Vec4f w_2 = { w.x * w.x, w.y * w.y, w.z * w.z };
	Vec4f w_4 = { w_2.x * w_2.x, w_2.y * w_2.y, w_2.z * w_2.z };
	Vec4f d_w = { d.z * w_4.z, d.z * w_4.z ,d.z * w_4.z };

	Vec4f five_two = { 52.0,52.0,52.0,52.0 };
	vector_time_matrix(s_i1, rot1, 1.0);
	return dotProduct4(d, five_two) * 0.5 + 0.5;
}

float Noise::simplex3d_fractal(Vec3f m)
{
	Vec3f m_eight = { 8 * m.x, 8 * m.y, 8 * m.z };
	return   0.5333333 * simplex3d(vector_time_matrix(m, rot1, 1.0))
		+ 0.2666667 * simplex3d(vector_time_matrix(m, rot2, 2.0))
		+ 0.1333333 * simplex3d(vector_time_matrix(m, rot3, 4.0))
		+ 0.666667 * simplex3d(m_eight);
}