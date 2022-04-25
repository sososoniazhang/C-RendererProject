/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"Matrix.h"
#include <iostream>
#include <algorithm>
#include <cmath>




void GzRender::findLightKs(GzColor color, float* E, float* Normal, GzCoord pointM) {
	GzColor sum = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < numlights; i++) {

		GzDepthMap* light = lights[i];

		
		if (light->GzFindZ(pointM) > 0.1) {
			continue;
		}
		


		GzMatrix Xiw;
		float N[3];
		float L[3];
		float dir[3];
		getUnitaryProtation(Xiw, m_camera.Xiw);
		for (int i = 0; i < 3; i++) {
			dir[i] = -(light->lightSource.lookat[i] - light->lightSource.position[i]);
		}
		project(Xiw, dir, L, 3);
		normalize(L);
		copyVector(Normal, N, 3);
		int lightContr = lightContribute(L, N, E);
		if (lightContr == -1) {
			flipNormal(N);
		}
		else if (lightContr == 0) {
			continue;
		}

		float R[3];
		findR(R, N, L);
		float RE = dot(R, E);

		if (RE <= 0) {
			RE = 0;
		}
		else if (RE > 1) {
			RE = 1;
		}
		float REs = pow(RE, spec);

		sum[0] += light->lightColor[0] * REs;
		sum[1] += light->lightColor[1] * REs;
		sum[2] += light->lightColor[2] * REs;
	}
	for (int i = 0; i < 3; i++) {
		if (sum[i] > 1) {
			sum[i] = 1.0;
		}
	}
	for (int i = 0; i < 3; i++) {
		(*(color + i)) = sum[i] * (Ks[i]);
		if (*(color + i) < 0) {
			*(color + i) = 0;
		}
		else if (*(color + i) > 1) {
			*(color + i) = 1;
		}
		else {
			continue;
		}
	}
	return;
}

void GzRender::findLightKd(GzColor color, float* E, float* Normal, GzCoord pointM) {
	GzColor sum = { 0.0, 0.0, 0.0 };
	GzColor kd_shade = { Kd[0], Kd[1], Kd[2] };

	for (int i = 0; i < numlights; i++) {
		GzDepthMap* light = lights[i];

		///*
		float shade = light->GzFindZ(pointM);
		if (shade > 0.1) {
			
						kd_shade[0] -= shade * Kd[0];
						kd_shade[1] -= shade * Kd[1];
						kd_shade[2] -= shade * Kd[2];
			

			continue;
		}
		//*/

		GzMatrix Xiw;
		float N[3];
		float L[3];
		float dir[3];
		getUnitaryProtation(Xiw, m_camera.Xiw);
		for (int i = 0; i < 3; i++) {
			dir[i] = -(light->lightSource.lookat[i] - light->lightSource.position[i]);
		}
		project(Xiw, dir, L, 3);
		normalize(L);
		copyVector(Normal, N, 3);
		int lightContr = lightContribute(L, N, E);
		if (lightContr == -1) {
			flipNormal(N);
		}
		else if (lightContr == 0) {
			continue;
		}


		float NLs = dot(N, L);

		if (NLs <= 0) {
			NLs = 0;
		}
		sum[0] += light->lightColor[0] * NLs;
		sum[1] += light->lightColor[1] * NLs;
		sum[2] += light->lightColor[2] * NLs;
	}
	for (int i = 0; i < 3; i++) {
		if (sum[i] > 1) {
			sum[i] = 1.0;
		}
	}


	for (int i = 0; i < 3; i++) {
		(*(color + i)) = sum[i] * (kd_shade[i]);
		if (*(color + i) < 0) {
			*(color + i) = 0;
		}
		else if (*(color + i) > 1) {
			*(color + i) = 1;
		}
		else {
			continue;
		}
	}
	return;
}



void interpolateVect(float* newVect, float* oldVect1, float* oldVect2, float* oldVect3, int x, int y) {
	for (int i = 0; i < 3; i++) {
		*(newVect + i) = generalLinearInt(*(oldVect1 + i), *(oldVect2 + i), *(oldVect3 + i), x, y);
	}
}

void GzRender::getLightSum(GzColor shade, GzCoord Normal, GzCoord pointM) {
	float specular[3] = { 0,0,0 };
	float defuse[3] = { 0,0,0 };
	GzCoord E;

	E[0] = 0.0;
	E[1] = 0.0;
	E[2] = -1.0;

	findLightKs(specular, E, Normal, pointM);
	findLightKd(defuse, E, Normal, pointM);
	for (int i = 0; i < 3; i++) {
		*(shade + i) = specular[i] + defuse[i] + Ka[i] * ambientlight.color[i];
		if (*(shade + i) > 1) {
			*(shade + i) = 1;
		}
	}
	return;
}
void GzRender::getColorFlat(GzColor shade, int i, int j, GzCoord pointM) {
	GzCoord norm;
	normalize(firstNorm);
	project(Xnorm[matlevel - 1], firstNorm, norm, 3);
	normalize(norm);

	getLightSum(shade, norm, pointM);

}
void GzRender::getColorGouraud(GzColor shade, int i, int j, GzCoord pointM) {
	GzCoord norm1, norm2, norm3;
	GzColor color1, color2, color3;
	project(Xnorm[matlevel - 1], myTriNormal[0], norm1, 3);
	project(Xnorm[matlevel - 1], myTriNormal[1], norm2, 3);
	project(Xnorm[matlevel - 1], myTriNormal[2], norm3, 3);
	normalize(norm1);
	normalize(norm2);
	normalize(norm3);

	getLightSum(color1, norm1, pointM);
	getLightSum(color2, norm2, pointM);
	getLightSum(color3, norm3, pointM);

	interpolateVect(shade, color1, color2, color3, i, j);

}
void GzRender::getColorPhong(GzColor shade, int i, int j, GzCoord pointM) {
	GzCoord norm1, norm2, norm3, normInt;
	GzCoord E = { 0, 0, -1 };
	project(Xnorm[matlevel - 1], myTriNormal[0], norm1, 3);
	project(Xnorm[matlevel - 1], myTriNormal[1], norm2, 3);
	project(Xnorm[matlevel - 1], myTriNormal[2], norm3, 3);
	normalize(norm1);
	normalize(norm2);
	normalize(norm3);
	interpolateVect(normInt, norm1, norm2, norm3, i, j);
	normalize(normInt);
	getLightSum(shade, normInt, pointM);

}


void GzRender::span(float* sx, float* sy, float* sz, float* ex, float* ey, float* ez, bool firstHalf) {
	if ((*sx == *ex)) {
		return;
	}
	float slopez = (*ez - *sz) / (*ex - *sx);
	float d_x = ceil(*sx) - *sx;
	int x, y;
	x = (int)ceil(*sx);
	y = (int)*sy;
	float z = *sz;
	while (x <= *ex) {
		z = z + slopez * d_x;
		x += 1;
		int index = y * yres + x;
		if (inBound(x, y, xres, yres) && z < pixelbuffer[index].z) {
			this->GzPut(x, y, ctoi(flatcolor[RED]), ctoi(flatcolor[GREEN]), ctoi(flatcolor[BLUE]), 0, (int)z);
		}

	}
}

void GzRender::advanceLine() {
	if (isFlat == 1) {
		float d_y = ceil(v1y) - v1y;
		moveDownStart(&startx, &starty, &startz, d_y, true);
		moveDownEnd(&endx, &endy, &endz, d_y, true);
		while (starty < v2y) {
			span(&startx, &starty, &startz, &endx, &endy, &endz, true);
			d_y = 1;
			moveDownStart(&startx, &starty, &startz, d_y, true);
			moveDownEnd(&endx, &endy, &endz, d_y, true);
		}
		return;
	}
	if (isFlat == -1) {
		startx = v1x;
		endx = v2x;
		startz = v1z;
		endz = v2z;
		float d_y = ceil(v1y) - v1y;
		moveDownStart(&startx, &starty, &startz, d_y, false);
		moveDownEnd(&endx, &endy, &endz, d_y, false);
		while (starty < v3y) {
			span(&startx, &starty, &startz, &endx, &endy, &endz, false);
			d_y = 1;
			moveDownStart(&startx, &starty, &startz, d_y, false);
			moveDownEnd(&endx, &endy, &endz, d_y, false);
		}
		return;
	}
	float d_y = ceil(v1y) - v1y;
	moveDownStart(&startx, &starty, &startz, d_y, true);
	moveDownEnd(&endx, &endy, &endz, d_y, true);
	while (starty <= v2y) {
		span(&startx, &starty, &startz, &endx, &endy, &endz, true);
		d_y = 1;
		moveDownStart(&startx, &starty, &startz, d_y, true);
		moveDownEnd(&endx, &endy, &endz, d_y, true);
	}
	d_y = v2y - starty;
	moveDownStart(&startx, &starty, &startz, d_y, true);
	moveDownEnd(&endx, &endy, &endz, d_y, true);

	// change to the second half
	d_y = ceil(v2y) - v2y;
	moveDownStart(&startx, &starty, &startz, d_y, false);
	moveDownEnd(&endx, &endy, &endz, d_y, false);
	while (starty <= v3y) {
		span(&startx, &starty, &startz, &endx, &endy, &endz, false);
		d_y = 1;
		moveDownStart(&startx, &starty, &startz, d_y, false);
		moveDownEnd(&endx, &endy, &endz, d_y, false);


	}
}


/***********************************************//***********************************************//***********************************************/
/***********************************************//***********************************************//***********************************************/
/***********************  start of HW 3 ************************/


int GzRender::GzRotXMat(float degree, GzMatrix mat)
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

int GzRender::GzRotYMat(float degree, GzMatrix mat)
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

int GzRender::GzRotZMat(float degree, GzMatrix mat)
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

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
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


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
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


GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */
	xres = xRes;
	yres = yRes;
	pixelbuffer = (GzPixel*)malloc(sizeof(GzPixel) * xres * yres);
	framebuffer = (char*)malloc(3 * sizeof(char) * xres * yres);
	numlights = 0;
	interp_mode = GZ_FLAT;

	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/
	matlevel = 0;
	m_camera.FOV = DEFAULT_FOV;
	m_camera.lookat[X] = 0.0;
	m_camera.lookat[Y] = 0.0;
	m_camera.lookat[Z] = 0.0;

	m_camera.position[X] = DEFAULT_IM_X;
	m_camera.position[Y] = -DEFAULT_IM_Y;
	m_camera.position[Z] = -DEFAULT_IM_Z;

	m_camera.worldup[X] = 0.0;
	m_camera.worldup[Y] = 1.0;
	m_camera.worldup[Z] = 0.0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				Xsp[i][j] = 1.0;
				m_camera.Xiw[i][j] = 1.0;
				m_camera.Xpi[i][j] = 1.0;
			}
			else {
				Xsp[i][j] = 0.0;
				m_camera.Xiw[i][j] = 0.0;
				m_camera.Xpi[i][j] = 0.0;
			}
		}
	}


}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	free(pixelbuffer);
	free(framebuffer);

}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */
	for (int i = 0; i < xres * yres; i++) {
		pixelbuffer[i].alpha = 4095;
		pixelbuffer[i].red = 3000;
		pixelbuffer[i].green = 3000;
		pixelbuffer[i].blue = 4000;
		pixelbuffer[i].z = INT_MAX;
	}

	return GZ_SUCCESS;
}


void GzRender::getXmp() {
	inverseMatrix(Ximage[matlevel - 1], Xmp);
}

int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/
	camDis = 1 / (tan((degToRadius(m_camera.FOV) / 2)));

	Xsp[0][0] = xres / 2;
	Xsp[0][3] = xres / 2;
	Xsp[1][1] = -xres / 2;
	Xsp[1][3] = yres / 2;
	Xsp[2][2] = MAXINT;
	GzPushMatrix(Xsp, 1);

	m_camera.Xpi[2][2] = 1 / camDis;
	m_camera.Xpi[3][2] = 1 / camDis;
	GzPushMatrix(m_camera.Xpi, 1);

	getXiw(&m_camera);
	GzPushMatrix(m_camera.Xiw);

	inverseMatrix(Ximage[matlevel - 1], Xwp);

	return GZ_SUCCESS;
}
int GzRender::GzPutShader(GzDepthMap* shadowMap) {
	lights[numlights] = shadowMap;
	numlights += 1;
	return GZ_SUCCESS;
}
int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/
	m_camera.FOV = camera.FOV;
	m_camera.lookat[X] = camera.lookat[X];
	m_camera.lookat[Y] = camera.lookat[Y];
	m_camera.lookat[Z] = camera.lookat[Z];

	m_camera.position[X] = camera.position[X];
	m_camera.position[Y] = camera.position[Y];
	m_camera.position[Z] = camera.position[Z];

	m_camera.worldup[X] = camera.worldup[X];
	m_camera.worldup[Y] = camera.worldup[Y];
	m_camera.worldup[Z] = camera.worldup[Z];

	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix, int stack)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	1:
		Ximage[mat], Xnormal[Mat], Xz[mat] Xz is world to perspective mat
	0:
		Ximage[mat], Xnormal[ID], Xz[ID]
	3:
		Ximage[mat], Xnormal[mat], Xz[ID]
	*/
	GzMatrix	pureRot =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};

	if (stack == 2) {
		getUnitaryProtation(pureRot, matrix);
	}
	if (matlevel == 0) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Ximage[0][i][j] = matrix[i][j];
				Xnorm[0][i][j] = pureRot[i][j];
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
				Xnorm[matlevel][i][j] = Xnorm[matlevel - 1][i][0] * pureRot[0][j] + Xnorm[matlevel - 1][i][1] * pureRot[1][j] + Xnorm[matlevel - 1][i][2] * pureRot[2][j] + Xnorm[matlevel - 1][i][3] * pureRot[3][j];

			}
		}
	}
	matlevel++;
	return GZ_SUCCESS;

}

int GzRender::GzPopMatrix(int stack)
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel == 0) {
		return GZ_FAILURE;
	}
	matlevel--;
	return GZ_SUCCESS;
}



int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
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


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
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


int GzRender::GzFlushDisplay2File(FILE* outfile)
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

int GzRender::GzFlushDisplay2FrameBuffer()
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


int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/
	for (int i = 0; i < numAttributes; i++) {
		if (*(nameList + i) == GZ_RGB_COLOR) {
			GzColor* value = (GzColor*)valueList[i];
			flatcolor[0] = (*value)[0];
			flatcolor[1] = (*value)[1];
			flatcolor[2] = (*value)[2];
		}
		else if (*(nameList + i) == GZ_AMBIENT_LIGHT) {
			GzLight* light = (GzLight*)valueList[i];
			ambientlight.color[0] = light->color[0];
			ambientlight.color[1] = light->color[1];
			ambientlight.color[2] = light->color[2];
			ambientlight.direction[0] = light->direction[0];
			ambientlight.direction[1] = light->direction[1];
			ambientlight.direction[2] = light->direction[2];
		}
		else if (*(nameList + i) == GZ_DIFFUSE_COEFFICIENT) {
			GzColor* value = (GzColor*)valueList[i];
			Kd[0] = (*value)[0];
			Kd[1] = (*value)[1];
			Kd[2] = (*value)[2];
		}
		else if (*(nameList + i) == GZ_AMBIENT_COEFFICIENT) {
			GzColor* value = (GzColor*)valueList[i];
			Ka[0] = (*value)[0];
			Ka[1] = (*value)[1];
			Ka[2] = (*value)[2];
		}
		else if (*(nameList + i) == GZ_SPECULAR_COEFFICIENT) {
			GzColor* value = (GzColor*)valueList[i];
			Ks[0] = (*value)[0];
			Ks[1] = (*value)[1];
			Ks[2] = (*value)[2];
		}
		else if (*(nameList + i) == GZ_DISTRIBUTION_COEFFICIENT) {
			spec = *(float*)valueList[i];
		}
		else if (*(nameList + i) == GZ_INTERPOLATE) {
			interp_mode = *(int*)valueList[i];
		}

	}
	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int	numParts, GzToken* nameList, GzPointer* valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/

	// get the coord and normal for each vertex

	for (int i = 0; i < numParts; i++) {
		if (nameList[i] == GZ_POSITION) {
			GzCoord* value = (GzCoord*)valueList[i];
			myTri[0][0] = (*value)[0];
			myTri[0][1] = (*value)[1];
			myTri[0][2] = (*value)[2];
			value++;
			myTri[1][0] = (*value)[0];
			myTri[1][1] = (*value)[1];
			myTri[1][2] = (*value)[2];
			value++;
			myTri[2][0] = (*value)[0];
			myTri[2][1] = (*value)[1];
			myTri[2][2] = (*value)[2];

		}
		if (nameList[i] == GZ_NORMAL) {
			GzCoord* value = (GzCoord*)valueList[i];
			myTriNormal[0][0] = (*value)[0];
			myTriNormal[0][1] = (*value)[1];
			myTriNormal[0][2] = (*value)[2];
			firstNorm[0] = (*value)[0];
			firstNorm[1] = (*value)[1];
			firstNorm[2] = (*value)[2];
			value++;
			myTriNormal[1][0] = (*value)[0];
			myTriNormal[1][1] = (*value)[1];
			myTriNormal[1][2] = (*value)[2];
			value++;
			myTriNormal[2][0] = (*value)[0];
			myTriNormal[2][1] = (*value)[1];
			myTriNormal[2][2] = (*value)[2];

		}
	}

	// transform to screen space
	float vector[4];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			vector[j] = myTri[i][0] * Ximage[matlevel - 1][j][0] + myTri[i][1] * Ximage[matlevel - 1][j][1] + myTri[i][2] * Ximage[matlevel - 1][j][2] + 1 * Ximage[matlevel - 1][j][3];
		}
		vector[3] = myTri[i][0] * Ximage[matlevel - 1][3][0] + myTri[i][1] * Ximage[matlevel - 1][3][1] + myTri[i][2] * Ximage[matlevel - 1][3][2] + 1 * Ximage[matlevel - 1][3][3];
		TriW[i] = vector[3];
		for (int j = 0; j < 3; j++) {
			myTri[i][j] = vector[j] / vector[3];
			if (j == 2 && myTri[i][j] < 0) {
				return GZ_FAILURE;
			}
		}

	}


	int* flag;
	if (sortVert(myTri, myTriNormal) == on_line_flag) {
		return GZ_SUCCESS;
	}
	// get the 3 vertex float number from sorted vertex
	getCoordTri(myTri);


	// find the bounding box
	int xStart = (int)floor(min(min((*myTri)[X], (*(myTri + 1))[X]), (*(myTri + 2))[X]));
	int xEnd = (int)ceil(max(max((*myTri)[X], (*(myTri + 1))[X]), (*(myTri + 2))[X]));
	int yStart = (int)floor(min(min((*myTri)[Y], (*(myTri + 1))[Y]), (*(myTri + 2))[Y]));
	int yEnd = (int)ceil(max(max((*myTri)[Y], (*(myTri + 1))[Y]), (*(myTri + 2))[Y]));
	// draw the lines in CW order
	drawLineCW(myTri);
	// for storing A, B, C, D for z interpolation
	float a, b, c, d;
	// find the a, b, c, d of the plane
	findABCD(&a, &b, &c, &d);

	for (int i = xStart; i <= xEnd; i++) {
		for (int j = yStart; j <= yEnd; j++) {
			// check in bound
			if (!inBound(i, j, xres, yres)) {
				continue;
			}
			// check if point inside triangle
			int flag1 = pointInLine(line1Params, i, j);
			int flag2 = pointInLine(line2Params, i, j);
			int flag3 = pointInLine(line3Params, i, j);
			if ((flag1 == right_flag && flag2 == right_flag && flag3 == right_flag) || (flag1 == on_line_flag) || (flag3 == on_line_flag && !isLRR) || (isFlat != 0 && checkflag == 1 && flag1 == right_flag) || (isFlat != 0 && checkflag == 3 && flag3 == right_flag) || (isFlat != 0 && checkflag == 2 && flag2 == right_flag)) {
				// interpolateZ and compare
				//float z = interpolateZ(i, j, &a, &b, &c, &d);

				float z = generalLinearInt(v1z, v2z, v3z, i, j);
				float z1 = interpolateZ(i, j, &a, &b, &c, &d);
				float w = generalLinearInt(TriW[0], TriW[1], TriW[2], i, j);
				if (z != z1) {
					return GZ_FAILURE;
				}
				//float z = interpolateZweight(i, j);
				int index = j * xres + i;
				if (z < pixelbuffer[index].z) {
					float pointP[4] = { i * w, j * w, z * w, w };
					float pointW[4] = { 0.0, 0.0, 0.0 , 0.0 };
					GzMatrix test;
					GzCoord model3dPoint;

					float pointM[4] = { 0.0, 0.0, 0.0 , 0.0 };
					project4d(Xmp, pointP, pointM);

					matMultiplication(Xmp, Ximage[5], test);
					for (int i = 0; i < 3; i++) {
						model3dPoint[i] = pointM[i] / pointM[3];
					}

					float shade[3];
					if (interp_mode == GZ_FLAT) {
						getColorFlat(shade, i, j, model3dPoint);
					}
					else if (interp_mode == GZ_COLOR) {
						getColorGouraud(shade, i, j, model3dPoint);
					}
					else {
						getColorPhong(shade, i, j, model3dPoint);
					}
					float r, g, b;
					r = ctoi(shade[RED]);
					g = ctoi(shade[GREEN]);
					b = ctoi(shade[BLUE]);
					this->GzPut(i, j, r, g, b, 0, (int)z);
				}
			}
			else {
				if (flag1 == -1 || flag2 == -1 || flag3 == -1) {
					return GZ_FAILURE;
				}

			}
		}

	}
	return GZ_SUCCESS;
}