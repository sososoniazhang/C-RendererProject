#pragma once
#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"Matrix.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#define PI (float) 3.14159265358979323846
bool LEE = true;


/***********************************************/
/* LEE Method */

int right_flag = 0;
int left_flag = 1;
int on_line_flag = -1;


// for storing the A, B, C parameters for the line of triangle sides
float line1Params[3];
float line2Params[3];
float line3Params[3];

GzCoord myTri[3];
GzCoord myTriNormal[3];
GzCoord firstNorm;
float TriW[3];

// storing the vertex coord
float v1z, v2z, v3z, v1x, v2x, v3x, v1y, v2y, v3y;



// this is to store the orientation of the triangle LLR or LRR
int isLRR = true;
// this is to store the orientation of the triangle with horizontal line
int isFlat = 0; // 0 when not flat 1 is top signle -1 is bottom single
int checkflag = 0;
/***********************************************/
/* Line Scale method */

// store the line y value range
float starty, endy;
float startx, endx;
float startz, endz;

float slope1x, slope1z;
float slope2x, slope2z;
float slope3x, slope3z;

/***********************************************/
/* distance from view plane to focal point */
float camDis;

/***********************************************/
/* general helper functions */
void matMultiplication(GzMatrix mat1, GzMatrix mat2, GzMatrix newmat) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			(*(newmat + i))[j] = mat1[i][0] * mat2[0][j] + mat1[i][1] * mat2[1][j] + mat1[i][2] * mat2[2][j] + mat1[i][3] * mat2[3][j];
		}
	}
}


bool inBound(int i, int j, int xres, int yres) {
	if (i < 0 || i >= xres || j < 0 || j >= yres) {
		return false;
	}
	return true;
}
bool gluInvertMatrixHelper(const float m[16], float invOut[16])
{
	float inv[16], det;
	int i;

	inv[0] = m[5] * m[10] * m[15] -
		m[5] * m[11] * m[14] -
		m[9] * m[6] * m[15] +
		m[9] * m[7] * m[14] +
		m[13] * m[6] * m[11] -
		m[13] * m[7] * m[10];

	inv[4] = -m[4] * m[10] * m[15] +
		m[4] * m[11] * m[14] +
		m[8] * m[6] * m[15] -
		m[8] * m[7] * m[14] -
		m[12] * m[6] * m[11] +
		m[12] * m[7] * m[10];

	inv[8] = m[4] * m[9] * m[15] -
		m[4] * m[11] * m[13] -
		m[8] * m[5] * m[15] +
		m[8] * m[7] * m[13] +
		m[12] * m[5] * m[11] -
		m[12] * m[7] * m[9];

	inv[12] = -m[4] * m[9] * m[14] +
		m[4] * m[10] * m[13] +
		m[8] * m[5] * m[14] -
		m[8] * m[6] * m[13] -
		m[12] * m[5] * m[10] +
		m[12] * m[6] * m[9];

	inv[1] = -m[1] * m[10] * m[15] +
		m[1] * m[11] * m[14] +
		m[9] * m[2] * m[15] -
		m[9] * m[3] * m[14] -
		m[13] * m[2] * m[11] +
		m[13] * m[3] * m[10];

	inv[5] = m[0] * m[10] * m[15] -
		m[0] * m[11] * m[14] -
		m[8] * m[2] * m[15] +
		m[8] * m[3] * m[14] +
		m[12] * m[2] * m[11] -
		m[12] * m[3] * m[10];

	inv[9] = -m[0] * m[9] * m[15] +
		m[0] * m[11] * m[13] +
		m[8] * m[1] * m[15] -
		m[8] * m[3] * m[13] -
		m[12] * m[1] * m[11] +
		m[12] * m[3] * m[9];

	inv[13] = m[0] * m[9] * m[14] -
		m[0] * m[10] * m[13] -
		m[8] * m[1] * m[14] +
		m[8] * m[2] * m[13] +
		m[12] * m[1] * m[10] -
		m[12] * m[2] * m[9];

	inv[2] = m[1] * m[6] * m[15] -
		m[1] * m[7] * m[14] -
		m[5] * m[2] * m[15] +
		m[5] * m[3] * m[14] +
		m[13] * m[2] * m[7] -
		m[13] * m[3] * m[6];

	inv[6] = -m[0] * m[6] * m[15] +
		m[0] * m[7] * m[14] +
		m[4] * m[2] * m[15] -
		m[4] * m[3] * m[14] -
		m[12] * m[2] * m[7] +
		m[12] * m[3] * m[6];

	inv[10] = m[0] * m[5] * m[15] -
		m[0] * m[7] * m[13] -
		m[4] * m[1] * m[15] +
		m[4] * m[3] * m[13] +
		m[12] * m[1] * m[7] -
		m[12] * m[3] * m[5];

	inv[14] = -m[0] * m[5] * m[14] +
		m[0] * m[6] * m[13] +
		m[4] * m[1] * m[14] -
		m[4] * m[2] * m[13] -
		m[12] * m[1] * m[6] +
		m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] +
		m[1] * m[7] * m[10] +
		m[5] * m[2] * m[11] -
		m[5] * m[3] * m[10] -
		m[9] * m[2] * m[7] +
		m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] -
		m[0] * m[7] * m[10] -
		m[4] * m[2] * m[11] +
		m[4] * m[3] * m[10] +
		m[8] * m[2] * m[7] -
		m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] +
		m[0] * m[7] * m[9] +
		m[4] * m[1] * m[11] -
		m[4] * m[3] * m[9] -
		m[8] * m[1] * m[7] +
		m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] -
		m[0] * m[6] * m[9] -
		m[4] * m[1] * m[10] +
		m[4] * m[2] * m[9] +
		m[8] * m[1] * m[6] -
		m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0)
		return false;

	det = 1.0 / det;

	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return true;
}

bool inverseMatrix(GzMatrix mat, GzMatrix inverse) {
	float old[16];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			old[i + 4 * j] = mat[i][j];
		}
	}
	float newmat[16];
	bool r = gluInvertMatrixHelper(old, newmat);
	if (!r) {
		return r;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			inverse[i][j] = newmat[i + 4 * j];
		}
	}
	return true;
}




void swapCord(GzCoord* vertexList, int a, int b) {
	float vx = (*(vertexList + a))[X];
	float vy = (*(vertexList + a))[Y];
	float vz = (*(vertexList + a))[Z];
	(*(vertexList + a))[X] = (*(vertexList + b))[X];
	(*(vertexList + a))[Y] = (*(vertexList + b))[Y];
	(*(vertexList + a))[Z] = (*(vertexList + b))[Z];
	(*(vertexList + b))[X] = vx;
	(*(vertexList + b))[Y] = vy;
	(*(vertexList + b))[Z] = vz;
}
void getCoordTri(GzCoord* vertexList) {
	v1z = (*vertexList)[Z];
	v2z = (*(vertexList + 1))[Z];
	v3z = (*(vertexList + 2))[Z];
	v1x = (*vertexList)[X];
	v2x = (*(vertexList + 1))[X];
	v3x = (*(vertexList + 2))[X];
	v1y = (*vertexList)[Y];
	v2y = (*(vertexList + 1))[Y];
	v3y = (*(vertexList + 2))[Y];
}

int sortVert(GzCoord* vertexList, GzCoord* normalList)
{
	float temp;
	int count = 0;
	for (int i = 0; i < 2; i++) {
		for (int p = 0; p < 2 - i; p++) {
			count += 1;
			if ((*(vertexList + p))[Y] > (*(vertexList + p + 1))[Y]) {
				swapCord(vertexList, p, p + 1);
				swapCord(normalList, p, p + 1);
				temp = TriW[p];
				TriW[p] = TriW[p + 1];
				TriW[p + 1] = temp;
			}
		}
	}

	float highx = (*vertexList)[X];
	float midx = (*(vertexList + 1))[X];
	float lowx = (*(vertexList + 2))[X];
	float highy = (*vertexList)[Y];
	float midy = (*(vertexList + 1))[Y];
	float lowy = (*(vertexList + 2))[Y];
	if ((highx == midx && midx == lowx) || (highy == midy && midy == lowy)) {
		// vertices in a line set flag to ignore later
		return on_line_flag;
	}
	if (highy == midy) {
		isFlat = -1;
		if (highx > midx) {

			swapCord(vertexList, 0, 1);
			swapCord(normalList, 0, 1);
			temp = TriW[0];
			TriW[0] = TriW[1];
			TriW[1] = temp;
			checkflag = 3;
		}
		else {
			checkflag = 1;
		}
		isLRR = true;
		return left_flag;
	}
	else if (lowy == midy) {
		isFlat = 1;
		if (lowx > midx) {
			swapCord(vertexList, 1, 2);
			swapCord(normalList, 1, 2);
			temp = TriW[2];
			TriW[2] = TriW[1];
			TriW[1] = temp;
			checkflag = 1;
		}
		else {
			checkflag = 2;
		}
		isLRR = true;
		return left_flag;

	}
	else {
		float m, b;
		if (highx == lowx) {
			if (midx > lowx) {
				isLRR = true;
				return left_flag;
			}
			if (midx < lowx) {
				isLRR = false;
				swapCord(vertexList, 1, 2);
				swapCord(normalList, 1, 2);
				temp = TriW[2];
				TriW[2] = TriW[1];
				TriW[1] = temp;
				return right_flag;
			}
		}
		m = (highy - lowy) / (highx - lowx);
		b = lowy - m * lowx;
		float intersection_x = (midy - b) / m;
		if (intersection_x < midx) {
			isLRR = true;
			return left_flag;
		}
		else if (intersection_x > midx) {
			if (LEE) {
				swapCord(vertexList, 1, 2);
				swapCord(normalList, 1, 2);
				temp = TriW[2];
				TriW[2] = TriW[1];
				TriW[1] = temp;
			}

			isLRR = false;
			return right_flag;
		}
		else {
			return on_line_flag;
		}
	}

}

/***********************************************/
/* LEE Method helper functions */

void drawLineCW(GzCoord* vertexList) {
	float dx, dy;
	// E1 (tail3->head1)
	dy = (v1y - v3y);
	dx = (v1x - v3x);
	line1Params[0] = dy;
	line1Params[1] = -dx;
	line1Params[2] = dx * v3y - dy * v3x;
	// E2 (1->2)
	dy = (v2y - v1y);
	dx = (v2x - v1x);
	line2Params[0] = dy;
	line2Params[1] = -dx;
	line2Params[2] = dx * v1y - dy * v1x;
	// E3 (2->3)
	dy = (v3y - v2y);
	dx = (v3x - v2x);
	line3Params[0] = dy;
	line3Params[1] = -dx;
	line3Params[2] = dx * v2y - dy * v2x;
}

int pointInLine(float lineParams[], int x, int y) {
	float check = lineParams[0] * x + lineParams[1] * y + lineParams[2];
	if (check > 0) {
		return left_flag;
	}
	else if (check <= 0) {
		return right_flag;
	}
	else {
		return on_line_flag;
	}
}

float determinate(float a, float b, float c, float d) {
	return a * d - c * b;
}


void crossProd(float* v1, float* v2, float* a, float* b, float* c) {
	float c1 = *(v1);
	float c2 = *(v1 + 1);
	float c3 = *(v1 + 2);
	float c4 = *(v2);
	float c5 = *(v2 + 1);
	float c6 = *(v2 + 2);

	*a = determinate(c2, c3, c5, c6);
	*b = -determinate(c1, c3, c4, c6);
	*c = determinate(c1, c2, c4, c5);

}

// to find and store a b c
void findABCD(float* a, float* b, float* c, float* d) {
	float vector1[3] = { v2x - v1x, v2y - v1y, v2z - v1z };
	float vector2[3] = { v3x - v2x, v3y - v2y, v3z - v2z };
	crossProd(vector1, vector2, a, b, c);
	*d = -(*a * v1x + *b * v1y + *c * v1z);
}


float interpolateZ(int x, int y, float* a, float* b, float* c, float* d) {
	return -(*a * x + *b * y + *d) / *c;
}

float interpolateZweight(int x, int y) {

	// https://codeplea.com/triangular-interpolation
	float w1, w2, w3;
	w1 = ((v2y - v3y) * (x - v3x) + (v3x - v2x) * (y - v3y)) / ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
	w2 = ((v3y - v1y) * (x - v3x) + (v1x - v3x) * (y - v3y)) / ((v2y - v3y) * (v1x - v3x) + (v3x - v2x) * (v1y - v3y));
	w3 = 1 - w1 - w2;
	float z = w1 * v1z + w1 * v2z + w3 * v3z;
	return z;

}

void fourDto3d(float* newV) {
	for (int i = 0; i < 3; i++) {
		*(newV + i) = *(newV + i) / *(newV + 3);
	}
}

void project(GzMatrix transX, GzCoord oldV, GzCoord newV, int num) {
	// project a 3D vector to new using the Xtransf (ASSUME NO Z SCALING!) num = 3 3d, num = 4 4d

	if (num == 3) {
		for (int j = 0; j < num; j++) {
			float t = 0;
			for (int i = 0; i < num; i++) {
				t += oldV[i] * transX[j][i];
			}
			*(newV + j) = t;
		}
		return;
	}
	float vector[4];
	for (int j = 0; j < 3; j++) {
		vector[j] = oldV[0] * transX[j][0] + oldV[1] * transX[j][1] + oldV[2] * transX[j][2] + 1 * transX[j][3];
	}
	vector[3] = oldV[0] * transX[3][0] + oldV[1] * transX[3][1] + oldV[2] * transX[3][2] + 1 * transX[3][3];

	for (int i = 0; i < 3; i++) {
		*(newV + i) = vector[i] / vector[3];
	}
	return;
}

void project4d(GzMatrix transX, float* oldV, float* newV) {
	// project a 3D vector to new using the Xtransf num = 3 3d, num = 4 4d
	for (int j = 0; j < 4; j++) {
		*(newV + j) = (*oldV) * transX[j][0] + *(oldV + 1) * transX[j][1] + *(oldV + 2) * transX[j][2] + *(oldV + 3) * transX[j][3];
	}

	return;
}



float generalLinearInt(float z1, float z2, float z3, int x, int y) {
	float vector1[3] = { v2x - v1x, v2y - v1y, z2 - z1 };
	float vector2[3] = { v3x - v2x, v3y - v2y, z3 - z2 };
	float a, b, c, d;
	crossProd(vector1, vector2, &a, &b, &c);
	d = -(a * v1x + b * v1y + c * z1);
	return -(a * x + b * y + d) / c;
}


/***********************************************/
/* Line scale Method helper functions */

void setSlope() {
	startx = v1x;
	starty = v1y;
	startz = v1z;
	endx = v1x;
	endy = v1y;
	endz = v1z;
	slope1x = (v2x - v1x) / (v2y - v1y);
	slope2x = (v3x - v2x) / (v3y - v2y);
	slope3x = (v3x - v1x) / (v3y - v1y);

	slope1z = (v2z - v1z) / (v2y - v1y);
	slope2z = (v3z - v2z) / (v3y - v2y);
	slope3z = (v3z - v1z) / (v3y - v1y);
}

void moveDownStart(float* x, float* y, float* z, float d_y, bool firstHalf) {
	if (!isLRR) {
		if (firstHalf) {
			startx = *x + slope1x * d_y;
			startz = *z + slope1z * d_y;
		}
		else {
			startx = *x + slope2x * d_y;
			startz = *z + slope2z * d_y;
		}
		starty = *y + d_y;
	}
	else {
		startx = *x + slope3x * d_y;
		startz = *z + slope3z * d_y;
		starty = *y + d_y;
	}
	*x = startx;
	*y = starty;
	*z = startz;

}

void moveDownEnd(float* x, float* y, float* z, float d_y, bool firstHalf) {

	if (isLRR) {
		if (firstHalf) {
			endx = *x + slope1x * d_y;
			endz = *z + slope1z * d_y;
		}
		else {
			endx = *x + slope2x * d_y;
			endz = *z + slope2z * d_y;
		}
		endy = *y + d_y;
	}
	else {
		endx = *x + slope3x * d_y;
		endz = *z + slope3z * d_y;
		endy = *y + d_y;
	}
	*x = endx;
	*y = endy;
	*z = endz;

}

/***********************  HW3 helper functions ************************/

float degToRadius(float degree) {
	return degree * PI / 180;
}

void normalize(GzCoord x) {
	float mag = sqrt((*x) * (*x) + (*(x + 1)) * (*(x + 1)) + (*(x + 2)) * (*(x + 2)));
	(*x) = (*x) / mag;
	(*(x + 1)) = (*(x + 1)) / mag;
	(*(x + 2)) = (*(x + 2)) / mag;
}

float dot(GzCoord x, GzCoord y) {
	return (*x) * (*y) + (*(x + 1)) * (*(y + 1)) + (*(x + 2)) * (*(y + 2));
}

void getXiwZ(GzCoord* z, GzCoord c, GzCoord l) {
	GzCoord cl;
	cl[0] = l[0] - c[0];
	cl[1] = *(l + 1) - *(c + 1);
	cl[2] = *(l + 2) - *(c + 2);
	normalize(cl);
	(*z)[0] = *(cl);
	(*z)[1] = *(cl + 1);
	(*z)[2] = *(cl + 2);
}

void getXiwY(GzCoord* y, GzCoord up, GzCoord z) {
	GzCoord upD;
	float dotP = dot(up, z);
	upD[0] = *(up)-(dot(up, z) * (*z));
	upD[1] = *(up + 1) - (dot(up, z) * (*(z + 1)));
	upD[2] = *(up + 2) - (dot(up, z) * (*(z + 2)));
	normalize(upD);
	(*y)[0] = *(upD);
	(*y)[1] = *(upD + 1);
	(*y)[2] = *(upD + 2);
}

void getXiwX(GzCoord* x, GzCoord y, GzCoord z) {
	float a, b, c;
	crossProd(y, z, &a, &b, &c);
	(*x)[0] = a;
	(*x)[1] = b;
	(*x)[2] = c;
}



void getXiw(GzCamera* cameraPtr) {
	GzCoord x, y, z;
	getXiwZ(&z, cameraPtr->position, cameraPtr->lookat);
	getXiwY(&y, cameraPtr->worldup, z);
	getXiwX(&x, y, z);


	(cameraPtr->Xiw)[0][0] = x[X];
	(cameraPtr->Xiw)[0][1] = x[Y];
	(cameraPtr->Xiw)[0][2] = x[Z];

	(cameraPtr->Xiw)[1][0] = y[X];
	(cameraPtr->Xiw)[1][1] = y[Y];
	(cameraPtr->Xiw)[1][2] = y[Z];

	(cameraPtr->Xiw)[2][0] = z[X];
	(cameraPtr->Xiw)[2][1] = z[Y];
	(cameraPtr->Xiw)[2][2] = z[Z];

	(cameraPtr->Xiw)[0][3] = -dot(x, cameraPtr->position);
	(cameraPtr->Xiw)[1][3] = -dot(y, cameraPtr->position);
	(cameraPtr->Xiw)[2][3] = -dot(z, cameraPtr->position);

}
/***********************  HW4 helper functions ************************/
/***********************  HW4 helper functions ************************/
/***********************  HW4 helper functions ************************/
/***********************  HW4 helper functions ************************/
/***********************  HW4 helper functions ************************/
/***********************  HW4 helper functions ************************/
void copyVector(float* fromV, float* toV, int len) {
	for (int i = 0; i < len; i++) {
		*(toV + i) = *(fromV + i);
	}
	return;
}


void findR(float* R, float* N, float* L) {
	float NL = dot(N, L);
	(*R) = 2 * NL * (*N) - (*L);
	(*(R + 1)) = 2 * NL * (*(N + 1)) - (*(L + 1));
	(*(R + 2)) = 2 * NL * (*(N + 2)) - (*(L + 2));
}


float getUnitaryProtation(GzMatrix uniRMat, GzMatrix mat) {

	float k = pow((mat[0][0] * mat[0][0] + mat[0][1] * mat[0][1] + mat[0][2] * mat[0][2]), 0.5);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			(*(uniRMat + i))[j] = mat[i][j] / k;
		}
	}
	return k;
}

int lightContribute(float* L, float* N, float* E) {
	float ln = dot(L, N);
	float en = dot(N, E);
	if (ln > 0 && en > 0) {
		return 1;
	}
	if (ln < 0 && en < 0) {
		return -1;
	}
	return 0;
}

void flipNormal(float* N) {
	// flip the normal vector
	for (int i = 0; i < 3; i++) {
		(*(N + i)) = -(*(N + i));
	}
	return;
}



void inverseRSTMatrix(GzMatrix mat, GzMatrix inverse) {
	GzMatrix rot =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};
	GzMatrix rot_inv =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};
	GzMatrix translation_inv =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};
	GzMatrix scale_inv =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};
	float s = getUnitaryProtation(rot, mat);
	for (int i = 0; i < 3; i++) {
		translation_inv[i][3] = -mat[i][3];
		scale_inv[i][i] = 1 / s;
		for (int j = 0; j < 3; j++) {
			rot_inv[i][j] = rot[j][i];
		}
	}
	GzMatrix st_inv =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};
	matMultiplication(scale_inv, translation_inv, st_inv);
	matMultiplication(rot_inv, st_inv, inverse);
	return;
}