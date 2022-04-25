#include "stdafx.h"
#include	"stdio.h"
#include "Gz.h"
#include "math.h"
#include "Zbuffer.h"
#include <stdlib.h>  
#include <stdio.h>  
#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;


/* LEE Method */

int right_flag_temp = 0;
int left_flag_temp = 1;
int on_line_flag_temp = -1;


// for storing the A, B, C parameters for the line of triangle sides
float line1Params_temp[3];
float line2Params_temp[3];
float line3Params_temp[3];


GzCoord myTri_temp[3];
GzCoord myTri_tempNormal[3];

// this is to store the orientation of the triangle LLR or LRR
int isLRR_temp = true;
int isFlat_temp = 0;





/***********************************************/
/* general helper functions */




void GzDepthMap::swapCord(GzCoord* vertexList, int a, int b) {
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
void GzDepthMap::getCoordTri(GzCoord* vertexList) {
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

int GzDepthMap::sortVert(GzCoord* vertexList, GzCoord* normalList)
{
	float temp;
	int count = 0;
	for (int i = 0; i < 2; i++) {
		for (int p = 0; p < 2 - i; p++) {
			count += 1;
			if ((*(vertexList + p))[Y] > (*(vertexList + p + 1))[Y]) {
				swapCord(vertexList, p, p + 1);
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
		return on_line_flag_temp;
	}
	if (highy == midy) {
		isFlat_temp = -1;
		if (highx > midx) {

			swapCord(vertexList, 0, 1);
		}
		isLRR_temp = true;
		return left_flag_temp;
	}
	else if (lowy == midy) {
		isFlat_temp = 1;
		if (lowx > midx) {
			swapCord(vertexList, 1, 2);
		}
		isLRR_temp = true;
		return left_flag_temp;

	}
	else {
		float m, b;
		if (highx == lowx) {
			if (midx > lowx) {
				isLRR_temp = true;
				return left_flag_temp;
			}
			if (midx < lowx) {
				isLRR_temp = false;
				swapCord(vertexList, 1, 2);
				return right_flag_temp;
			}
		}
		m = (highy - lowy) / (highx - lowx);
		b = lowy - m * lowx;
		float intersection_x = (midy - b) / m;
		if (intersection_x < midx) {
			isLRR_temp = true;
			return left_flag_temp;
		}
		else if (intersection_x > midx) {
			swapCord(vertexList, 1, 2);

			isLRR_temp = false;
			return right_flag_temp;
		}
		else {
			return on_line_flag_temp;
		}
	}

}


/***********************************************/
/* LEE Method helper functions */

void GzDepthMap::drawLineCW(GzCoord* vertexList) {
	float dx, dy;
	// E1 (tail3->head1)
	dy = (v1y - v3y);
	dx = (v1x - v3x);
	line1Params_temp[0] = dy;
	line1Params_temp[1] = -dx;
	line1Params_temp[2] = dx * v3y - dy * v3x;
	// E2 (1->2)
	dy = (v2y - v1y);
	dx = (v2x - v1x);
	line2Params_temp[0] = dy;
	line2Params_temp[1] = -dx;
	line2Params_temp[2] = dx * v1y - dy * v1x;
	// E3 (2->3)
	dy = (v3y - v2y);
	dx = (v3x - v2x);
	line3Params_temp[0] = dy;
	line3Params_temp[1] = -dx;
	line3Params_temp[2] = dx * v2y - dy * v2x;
}

int GzDepthMap::pointInLine(float lineParams[], int x, int y) {
	float check = lineParams[0] * x + lineParams[1] * y + lineParams[2];
	if (check > 0) {
		return left_flag_temp;
	}
	else if (check < 0) {
		return right_flag_temp;
	}
	else {
		return on_line_flag_temp;
	}
}


void GzDepthMap::crossProd(float* v1, float* v2, float* a, float* b, float* c) {
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
void GzDepthMap::findABCD(float* a, float* b, float* c, float* d) {
	float vector1[3] = { v2x - v1x, v2y - v1y, v2z - v1z };
	float vector2[3] = { v3x - v2x, v3y - v2y, v3z - v2z };
	crossProd(vector1, vector2, a, b, c);
	*d = -(*a * v1x + *b * v1y + *c * v1z);
}




void GzDepthMap::fourDto3d(float* newV) {
	for (int i = 0; i < 3; i++) {
		*(newV + i) = *(newV + i) / *(newV + 3);
	}
}

void GzDepthMap::project(GzMatrix transX, GzCoord oldV, GzCoord newV, int num) {
	// project a 3D vector to new using the Xtransf (ASSUME NO Z SCALING!) num = 3 3d, num = 4 4d
	for (int j = 0; j < num; j++) {
		float t = 0;
		for (int i = 0; i < num; i++) {
			t += oldV[i] * transX[j][i];
		}
		*(newV + j) = t;
	}
	return;
}



float GzDepthMap::generalLinearInt(float z1, float z2, float z3, int x, int y) {
	float vector1[3] = { v2x - v1x, v2y - v1y, z2 - z1 };
	float vector2[3] = { v3x - v2x, v3y - v2y, z3 - z2 };
	float a, b, c, d;
	crossProd(vector1, vector2, &a, &b, &c);
	d = -(a * v1x + b * v1y + c * z1);
	return -(a * x + b * y + d) / c;
}



void GzDepthMap::normalize(GzCoord x) {
	float mag = sqrt((*x) * (*x) + (*(x + 1)) * (*(x + 1)) + (*(x + 2)) * (*(x + 2)));
	(*x) = (*x) / mag;
	(*(x + 1)) = (*(x + 1)) / mag;
	(*(x + 2)) = (*(x + 2)) / mag;
}


void GzDepthMap::getXiwZ(GzCoord* z, GzCoord c, GzCoord l) {
	GzCoord cl;
	cl[0] = l[0] - c[0];
	cl[1] = *(l + 1) - *(c + 1);
	cl[2] = *(l + 2) - *(c + 2);
	normalize(cl);
	(*z)[0] = *(cl);
	(*z)[1] = *(cl + 1);
	(*z)[2] = *(cl + 2);
}

void GzDepthMap::getXiwY(GzCoord* y, GzCoord up, GzCoord z) {
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

void GzDepthMap::getXiwX(GzCoord* x, GzCoord y, GzCoord z) {
	float a, b, c;
	crossProd(y, z, &a, &b, &c);
	(*x)[0] = a;
	(*x)[1] = b;
	(*x)[2] = c;
}



void GzDepthMap::getXiw(GzCamera* cameraPtr) {
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

float GzDepthMap::interpolateDepth2d(float x, float y, int layer) {
	int xlow = floor(x);
	int ylow = floor(y);
	float result = -1.0;
	if (inBound(xlow, ylow, xres, yres) && inBound(xlow + 1, ylow + 1, xres, yres)) {
		float s = x - xlow;
		float t = y - ylow;
		float a, b, c, d;
		if (layer == 0) {
			a = zbuffer[ARRAY(xlow, ylow)];
			b = zbuffer[ARRAY(xlow + 1, ylow)];
			c = zbuffer[ARRAY(xlow + 1, ylow + 1)];
			d = zbuffer[ARRAY(xlow, ylow + 1)];
		}
		else {
			a = zbuffer[ARRAY(xlow, ylow)];
			b = zbuffer[ARRAY(xlow + 1, ylow)];
			c = zbuffer[ARRAY(xlow + 1, ylow + 1)];
			d = zbuffer[ARRAY(xlow, ylow + 1)];
		}
		result = s * t * d + (1 - s) * t * c + (1 - t) * s * b + (1 - s) * (1 - t) * a;
	}
	return result;
}
float GzDepthMap::interpolateDepth3d(float x, float y, float z) {
	int zfloor = floor(z / step);
	float result = -1.0;
	if (layer == 1) {
		return interpolateDepth2d(x, y);
	}
	if (z < zlow) {
		return 0.0;
	}
	if (z > zhigh) {
		return 1.0;
	}
	float zsmall = interpolateDepth2d(x, y, zfloor);
	float zlarge = interpolateDepth2d(x, y, zfloor + 1);
	float ratio = z - zfloor * step;
	result = (ratio / step) * zlarge + (1 - (ratio / step)) * zsmall;
	return result;
}

void GzDepthMap::putColor(GzColor color) {
	for (int i = 0; i < 3; i++) {
		lightColor[i] = color[i];
	}
}

int GzDepthMap::transWS(vertex* position) {
	// using Ximage[2] since we only need the transformation from world to perspective space (stack-> Xsp Xpi Xiw)
	float vector[4];
	for (int j = 0; j < 3; j++) {
		vector[j] = position->x * Ximage[matlevel - 1][j][0] + position->y * Ximage[matlevel - 1][j][1] + position->z * Ximage[matlevel - 1][j][2] + 1 * Ximage[matlevel - 1][j][3];
	}
	vector[3] = position->x * Ximage[matlevel - 1][3][0] + position->y * Ximage[matlevel - 1][3][1] + position->z * Ximage[matlevel - 1][3][2] + 1 * Ximage[matlevel - 1][3][3];

	for (int j = 0; j < 3; j++) {
		vector[j] = vector[j] / vector[3];
	}
	position->x = vector[0];
	position->y = vector[1];
	position->z = vector[2];
	return GZ_SUCCESS;
}

int GzDepthMap::readTriangle(vertex a, vertex b, vertex c) {
	myTri_temp[0][0] = a.x;
	myTri_temp[0][1] = a.y;
	myTri_temp[0][2] = a.z;

	myTri_temp[1][0] = b.x;
	myTri_temp[1][1] = b.y;
	myTri_temp[1][2] = b.z;

	myTri_temp[2][0] = c.x;
	myTri_temp[2][1] = c.y;
	myTri_temp[2][2] = c.z;


	return GZ_SUCCESS;
}

int GzDepthMap::drawMap() {
	// find the bounding box
	int xStart = (int)floor(min(min((*myTri_temp)[X], (*(myTri_temp + 1))[X]), (*(myTri_temp + 2))[X]));
	int xEnd = (int)ceil(max(max((*myTri_temp)[X], (*(myTri_temp + 1))[X]), (*(myTri_temp + 2))[X]));
	int yStart = (int)floor(min(min((*myTri_temp)[Y], (*(myTri_temp + 1))[Y]), (*(myTri_temp + 2))[Y]));
	int yEnd = (int)ceil(max(max((*myTri_temp)[Y], (*(myTri_temp + 1))[Y]), (*(myTri_temp + 2))[Y]));

	for (int i = xStart; i <= xEnd; i++) {
		for (int j = yStart; j <= yEnd; j++) {
			// check in bound
			if (!inBound(i, j, xres, yres)) {
				continue;
			}
			// check if point inside triangle
			int flag1 = pointInLine(line1Params_temp, i, j);
			int flag2 = pointInLine(line2Params_temp, i, j);
			int flag3 = pointInLine(line3Params_temp, i, j);
			if ((flag1 == right_flag_temp && flag2 == right_flag_temp && flag3 == right_flag_temp) || (flag1 == on_line_flag_temp) || (flag3 == on_line_flag_temp && !isLRR_temp)) {
				// interpolateZ and compare
				float z = interpolateZ(i, j, &A, &B, &C, &D);
				z += 10000000;
				if (zbuffer[ARRAY(i, j)] > z) {
					zbuffer[ARRAY(i, j)] = z;
				}
			}
		}
	}
	return GZ_SUCCESS;
}



GzDepthMap::GzDepthMap(int zmax, bool d) {
	threeD = d;
	xres = 3000;
	yres = 3000;
	maxZ = zmax;
	pixelbuffer = (GzPixel*)malloc(sizeof(GzPixel) * xres * yres);
	framebuffer = (char*)malloc(3 * sizeof(char) * xres * yres);
	zbuffer = (float*)malloc(sizeof(float) * xres * yres);
	for (int i = 0; i < xres * yres; i++) {
		zbuffer[i] = maxZ;
	}
	// SETUP default light source negative x from default camera?
	matlevel = 0;
	lightSource.FOV = DEFAULT_FOV;
	lightSource.lookat[X] = 0.0;
	lightSource.lookat[Y] = 0.0;
	lightSource.lookat[Z] = 0.0;

	lightSource.position[X] = DEFAULT_IM_X;
	lightSource.position[Y] = -DEFAULT_IM_Y;
	lightSource.position[Z] = -DEFAULT_IM_Z;

	lightSource.worldup[X] = 0.0;
	lightSource.worldup[Y] = 1.0;
	lightSource.worldup[Z] = 0.0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				Xsp[i][j] = 1.0;
				lightSource.Xiw[i][j] = 1.0;
				lightSource.Xpi[i][j] = 1.0;
			}
			else {
				Xsp[i][j] = 0.0;
				lightSource.Xiw[i][j] = 0.0;
				lightSource.Xpi[i][j] = 0.0;
			}
		}
	}
	for (int i = 0; i < xres * yres; i++) {
		pixelbuffer[i].alpha = 0;
		pixelbuffer[i].red = 4095;
		pixelbuffer[i].green = 4095;
		pixelbuffer[i].blue = 4095;
		pixelbuffer[i].z = maxZ;
	}
}
GzDepthMap::~GzDepthMap() {
	free(zbuffer);
	zbuffer = nullptr;
}


int GzDepthMap::GzPushMatrix(GzMatrix matrix)
{
	/*
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	if (matlevel == 0) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Ximage[0][i][j] = matrix[i][j];
			}
		}
	}
	else if (matlevel == MATLEVELS) {
		return GZ_FAILURE;
	}
	else {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Ximage[matlevel][i][j] = Ximage[matlevel - 1][i][0] * matrix[0][j] + Ximage[matlevel - 1][i][1] * matrix[1][j] + Ximage[matlevel - 1][i][2] * matrix[2][j] + Ximage[matlevel - 1][i][3] * matrix[3][j];
			}
		}
	}
	matlevel++;
	return GZ_SUCCESS;
}

int GzDepthMap::GzPopMatrix()
{
	/*
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel == 0) {
		return GZ_FAILURE;
	}
	matlevel--;
	return GZ_SUCCESS;
}


int GzDepthMap::GzPutLightSource(GzCamera lightsource)
{
	/*
	/*- overwrite renderer light source with new light definition
	*/
	lightSource.FOV = lightsource.FOV;
	lightSource.lookat[X] = lightsource.lookat[X];
	lightSource.lookat[Y] = lightsource.lookat[Y];
	lightSource.lookat[Z] = lightsource.lookat[Z];

	lightSource.position[X] = lightsource.position[X];
	lightSource.position[Y] = lightsource.position[Y];
	lightSource.position[Z] = lightsource.position[Z];

	lightSource.worldup[X] = lightsource.worldup[X];
	lightSource.worldup[Y] = lightsource.worldup[Y];
	lightSource.worldup[Z] = lightsource.worldup[Z];

	return GZ_SUCCESS;
}

int GzDepthMap::startUp() {

	camDis = 1 / (tan((degToRadius(lightSource.FOV) / 2)));

	Xsp[0][0] = xres / 2;
	Xsp[0][3] = xres / 2;
	Xsp[1][1] = -(float)(yres / 2);
	Xsp[1][3] = yres / 2;
	Xsp[2][2] = maxZ;
	GzPushMatrix(Xsp);

	lightSource.Xpi[2][2] = 1 / camDis;
	lightSource.Xpi[3][2] = 1 / camDis;
	GzPushMatrix(lightSource.Xpi);

	getXiw(&lightSource);
	GzPushMatrix(lightSource.Xiw);
	return GZ_SUCCESS;
}


// generate depth map from triangle
int GzDepthMap::GzCreateDepthMapFromSTL(StlObject* obj) {
	objptr = obj;
	for (int i = 0; i < obj->nTri; i++) {

		//get three vertex
		vertex a, b, c, n;
		obj->getTri(i, &a, &b, &c, &n);

		// read triangle
		readTriangle(a, b, c);
		// transform to screen space
		float vector[4];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				vector[j] = myTri_temp[i][0] * Ximage[matlevel - 1][j][0] + myTri_temp[i][1] * Ximage[matlevel - 1][j][1] + myTri_temp[i][2] * Ximage[matlevel - 1][j][2] + 1 * Ximage[matlevel - 1][j][3];
			}
			vector[3] = myTri_temp[i][0] * Ximage[matlevel - 1][3][0] + myTri_temp[i][1] * Ximage[matlevel - 1][3][1] + myTri_temp[i][2] * Ximage[matlevel - 1][3][2] + 1 * Ximage[matlevel - 1][3][3];

			for (int j = 0; j < 3; j++) {
				myTri_temp[i][j] = vector[j] / vector[3];
				if (j == 2 && myTri_temp[i][j] < 0) {
					return GZ_FAILURE;
				}
			}

		}
		int* flag;
		if (sortVert(myTri_temp, myTri_tempNormal) == on_line_flag_temp) {
			return GZ_SUCCESS;
		}
		// get the 3 vertex float number from sorted vertex
		getCoordTri(myTri_temp);

		// draw the lines in CW order
		drawLineCW(myTri_temp);

		// find the a, b, c, d of the plane
		findABCD(&A, &B, &C, &D);

		drawMap();
	}
	return GZ_SUCCESS;
}

int GzDepthMap::GzCreateDepthMapFromTri(int	numParts, GzToken* nameList, GzPointer* valueList) {
	for (int i = 0; i < numParts; i++) {
		if (nameList[i] == GZ_POSITION) {
			GzCoord* value = (GzCoord*)valueList[i];
			myTri_temp[0][0] = (*value)[0];
			myTri_temp[0][1] = (*value)[1];
			myTri_temp[0][2] = (*value)[2];
			value++;
			myTri_temp[1][0] = (*value)[0];
			myTri_temp[1][1] = (*value)[1];
			myTri_temp[1][2] = (*value)[2];
			value++;
			myTri_temp[2][0] = (*value)[0];
			myTri_temp[2][1] = (*value)[1];
			myTri_temp[2][2] = (*value)[2];

		}
	}
	// transform three vertex to screen space
	float vector[4];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			vector[j] = myTri_temp[i][0] * Ximage[matlevel - 1][j][0] + myTri_temp[i][1] * Ximage[matlevel - 1][j][1] + myTri_temp[i][2] * Ximage[matlevel - 1][j][2] + 1 * Ximage[matlevel - 1][j][3];
		}
		vector[3] = myTri_temp[i][0] * Ximage[matlevel - 1][3][0] + myTri_temp[i][1] * Ximage[matlevel - 1][3][1] + myTri_temp[i][2] * Ximage[matlevel - 1][3][2] + 1 * Ximage[matlevel - 1][3][3];

		for (int j = 0; j < 3; j++) {
			myTri_temp[i][j] = vector[j] / vector[3];
			if (j == 2 && myTri_temp[i][j] < 0) {
				return GZ_FAILURE;
			}
		}

	}
	int* flag;
	if (sortVert(myTri_temp, myTri_tempNormal) == on_line_flag_temp) {
		return GZ_SUCCESS;
	}
	// get the 3 vertex float number from sorted vertex
	getCoordTri(myTri_temp);

	// draw the lines in CW order
	drawLineCW(myTri_temp);

	// find the a, b, c, d of the plane
	findABCD(&A, &B, &C, &D);
	drawMap();

	return GZ_SUCCESS;
}


float GzDepthMap::GzFindZ(GzCoord position) {
	vertex point;
	point.x = position[0];
	point.y = position[1];
	point.z = position[2];
	transWS(&point);

	if (threeD) {
		return interpolateDepth3d(point.x, point.y, point.z);
	}
	float depth = interpolateDepth2d(point.x, point.y);
	if (depth < 0) {
		return -1;
	}
	if (depth < point.z) {
		return max(point.z - depth, (float)100) / 100; // shade
	}
	return 0.0; // no shade

}

int GzDepthMap::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	fprintf(outfile, "P6 %d %d %d\r", xres, yres, 255);
	for (int i = 0; i < xres * yres; i++) {
		int r = pixelbuffer[i].red >> 4;
		int g = pixelbuffer[i].green >> 4;
		int b = pixelbuffer[i].blue >> 4;
		fprintf(outfile, "%c%c%c", r, g, b);
	}


	return GZ_SUCCESS;
}

int GzDepthMap::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	int index = 0;
	for (int i = 0; i < xres * yres; i++) {
		int r = pixelbuffer[i].red >> 4;
		int g = pixelbuffer[i].green >> 4;
		int b = pixelbuffer[i].blue >> 4;
		framebuffer[index++] = b;
		framebuffer[index++] = g;
		framebuffer[index++] = r;
	}

	return GZ_SUCCESS;
}


int GzDepthMap::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */
	if (i < 0 || i >= xres || j < 0 || j >= yres) {
		return GZ_FAILURE;
	}
	int index = i + j * xres;

	pixelbuffer[index].alpha = (a < 0) ? 0 : a;
	pixelbuffer[index].alpha = (pixelbuffer[index].alpha > 4095) ? 4095 : pixelbuffer[index].alpha;

	pixelbuffer[index].red = (r < 0) ? 0 : r;
	pixelbuffer[index].red = (pixelbuffer[index].red > 4095) ? 4095 : pixelbuffer[index].red;
	pixelbuffer[index].green = (g < 0) ? 0 : g;
	pixelbuffer[index].green = (pixelbuffer[index].green > 4095) ? 4095 : pixelbuffer[index].green;
	pixelbuffer[index].blue = (b < 0) ? 0 : b;
	pixelbuffer[index].blue = (pixelbuffer[index].blue > 4095) ? 4095 : pixelbuffer[index].blue;
	pixelbuffer[index].z = (z < 0) ? 0 : z;
	pixelbuffer[index].z = (pixelbuffer[index].z > INT_MAX) ? INT_MAX : pixelbuffer[index].z;

	return GZ_SUCCESS;
}


int GzDepthMap::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */

	if (i < 0 || i >= xres || j < 0 || j >= yres) {
		return GZ_FAILURE;
	}
	int index = i + j * xres;
	*a = pixelbuffer[index].alpha;
	*r = pixelbuffer[index].red;
	*g = pixelbuffer[index].green;
	*b = pixelbuffer[index].blue;
	*z = pixelbuffer[index].z;
	return GZ_SUCCESS;
}


int GzDepthMap::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	float rad = degToRadius(degree);
	(*mat)[0] = 1;
	(*(mat + 1))[1] = cos(rad);
	(*(mat + 1))[2] = -sin(rad);
	(*(mat + 2))[1] = sin(rad);
	(*(mat + 2))[2] = cos(rad);

	return GZ_SUCCESS;
}

int GzDepthMap::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	float rad = degToRadius(degree);
	(*(mat))[0] = cos(rad);
	(*(mat))[2] = sin(rad);
	(*(mat + 2))[0] = -sin(rad);
	(*(mat + 2))[2] = cos(rad);

	return GZ_SUCCESS;
}

int GzDepthMap::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	float rad = degToRadius(degree);
	(*mat)[0] = cos(rad);
	(*mat)[1] = -sin(rad);
	(*(mat + 1))[0] = sin(rad);
	(*(mat + 1))[1] = cos(rad);

	return GZ_SUCCESS;
}

int GzDepthMap::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/
	(*mat)[3] = *translate;
	(*(mat + 1))[3] = *(translate + 1);
	(*(mat + 2))[3] = *(translate + 2);

	return GZ_SUCCESS;
}


int GzDepthMap::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	(*mat)[0] = *scale;
	(*(mat + 1))[1] = *(scale + 1);
	(*(mat + 2))[2] = *(scale + 2);

	return GZ_SUCCESS;
}

void GzDepthMap::getShadowMap() {
	float minz = maxZ;
	float maxz = -100000;
	for (int i = 0; i < xres * yres; i++) {
		if (zbuffer[i] == maxZ) {
			continue;
		}
		else {
			minz = min(minz, zbuffer[i]);
			maxz = max(maxz, zbuffer[i]);
		}
	}
	float range = maxz - minz;
	for (int i = 0; i < xres; i++) {
		for (int j = 0; j < yres; j++) {
			float z = zbuffer[ARRAY(i, j)];

			if (z == maxZ) {
				continue;
			}
			if (z < 0) {
				int i;
			}
			else {
				float r = min(4095, ctoi((z - minz) / range));
				this->GzPut(i, j, r, r, r, 4095, (int)z);
				float test = framebuffer[ARRAY(i, j)];
			}
		}
	}
}
