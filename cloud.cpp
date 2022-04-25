#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <math.h>
#include "cloud.h"
#include <algorithm>
#include "noise.h"
#define PI (float) 3.14159265358979323846

#define coef_C  0.01

using namespace std;

// Create the cloud object center and default value for camera
Cloud::Cloud(float centerX, float centerY, float centerZ)
{
	//Define the center of the cloud
	center_x = centerX;
	center_y = centerY;
	center_z = centerZ;

	/*Set up the default camera*/
	c_camera.lookat[X] = 0;
	c_camera.lookat[Y] = 0;
	c_camera.lookat[Z] = 0;

	c_camera.position[X] = DEFAULT_IM_X;
	c_camera.position[Y] = DEFAULT_IM_Y;
	c_camera.position[Z] = DEFAULT_IM_Z;

	c_camera.worldup[X] = 0;
	c_camera.worldup[Y] = 1;
	c_camera.worldup[Z] = 0;

	c_camera.FOV = DEFAULT_FOV;

}

Cloud::~Cloud() {};
/*
int Cloud::CloudDefault()
{	//Default sphere cloud
	sphereCloud.radius = 2.0;
	//Default cylinder cloud
	cylinderCloud.radius = 2.0;
	cylinderCloud.height = 1.5;

	//Default block cloud
	blockCloud.length = 2.5;
	blockCloud.width = 1.0;
	blockCloud.height = 1.5;

	//Default cube cloud
	cubeCloud.length = 2.0;

	return GZ_SUCCESS;
}
*/
int Cloud::putSphere(float radius)
{
	if (object_type == 0) {
		sphereCloud.radius = radius;
		object_type = SPHERE;

	}
	return GZ_SUCCESS;
}

int Cloud::putCylinder(float radius, float height)
{
	if (object_type == 0) {
		cylinderCloud.radius = radius;
		cylinderCloud.height = height;
		object_type = CYLINDER;

	}
	return GZ_SUCCESS;
}

int Cloud::putBlock(float length, float width, float height)
{
	if (object_type == 0) {
		blockCloud.length = length;
		blockCloud.width = width;
		blockCloud.height = height;
		object_type = BLOCK;
	}
	return GZ_SUCCESS;
}

int Cloud::putCube(float length)
{
	if (object_type == 0) {
		cubeCloud.length = length;
		object_type = CUBE;
	}
	return GZ_SUCCESS;
}


float Cloud::isCloud(GzCoord point)
{
	float rate = 0.0;
	if (object_type == SPHERE)
	{
		float distance = 0.0; // Distance from given point to the center of the cloud

		distance = sqrtf(powf(point[0] - center_x, 2.0) + powf(point[1] - center_y, 2.0) + powf(point[2] - center_z, 2.0));

		if (distance > sphereCloud.radius) {
			return rate;
		}
		else
		{
			Vec3f p = { point[0] * coef_C * 8.0 + 8.0,point[1] * coef_C * 8.0 + 8.0, point[2] * coef_C * 8.0 + 8.0 };
			rate = noise->simplex3d_fractal(p);
			rate = max(rate, 0.0);
			rate = min(rate, 1.0);
			return rate;
			//return 1.0;
		}
	}


	else if (object_type == BLOCK)
	{
		float distance_x = abs(center_x - point[0]);
		float distance_y = abs(center_y - point[1]);
		float distance_z = abs(center_z - point[2]);

		if (distance_x > (blockCloud.length / 2.0) || distance_y > (blockCloud.width / 2.0) || distance_z > (blockCloud.height / 2.0))
		{
			return rate;
		}
		else
		{
			Vec3f p = { point[0] * coef_C * 8.0 + 8.0,point[1] * coef_C * 8.0 + 8.0, point[2] * coef_C * 8.0 + 8.0 };
			rate = noise->simplex3d_fractal(p);
			rate = max(rate, 0.0);
			rate = min(rate, 1.0);
			return rate;
		}
	}

	else if (object_type == CUBE)
	{
		float distance_x = abs(center_x - point[0]);
		float distance_y = abs(center_y - point[1]);
		float distance_z = abs(center_z - point[2]);

		if (distance_x > (cubeCloud.length / 2.0) || distance_y > (cubeCloud.length / 2.0) || distance_z > (cubeCloud.length / 2.0))
		{
			return rate;
		}
		else
		{
			Vec3f p = { point[0] * coef_C * 8.0 + 8.0,point[1] * coef_C * 8.0 + 8.0, point[2] * coef_C * 8.0 + 8.0 };
			rate = noise->simplex3d_fractal(p);
			rate = max(rate, 0.0);
			rate = min(rate, 1.0);
			return rate;
			return rate;
		}
	}

	else if (object_type == CYLINDER)
	{
		float distance_x = abs(center_x - point[0]);
		float yz_center = sqrtf(powf(center_y - point[1], 2.0) + powf(center_z - point[2], 2.0));

		if (distance_x > (cylinderCloud.height / 2.0) || yz_center > cylinderCloud.radius)
		{
			return rate;
		}
		else
		{
			Vec3f p = { point[0] * coef_C * 8.0 + 8.0,point[1] * coef_C * 8.0 + 8.0, point[2] * coef_C * 8.0 + 8.0 };
			rate = noise->simplex3d_fractal(p);
			rate = max(rate, 0.0);
			rate = min(rate, 1.0);
			return rate;
		}
	}

}
//Transform from model to screen space
void transform(GzMatrix mat1, float* ver, float* result)
{
	float ver_screen[4][1];

	for (int i = 0; i < 4; i++) {
		ver_screen[i][0] = 0.0;
		for (int j = 0; j < 4; j++) {
			ver_screen[i][0] += (mat1[i][j] * ver[j]);
		}
	}
	result[0] = ver_screen[0][0] / ver_screen[3][0];
	result[1] = ver_screen[1][0] / ver_screen[3][0];
	result[2] = ver_screen[2][0] / ver_screen[3][0];

}

float* calculateZInterval(GzMatrix mat, boundary box)
{
	float ver1[4] = { box.v1.x, box.v1.y, box.v1.z, 1.0 };
	float ver2[4] = { box.v2.x, box.v2.y, box.v2.z, 1.0 };
	float ver3[4] = { box.v3.x, box.v3.y, box.v3.z, 1.0 };
	float ver4[4] = { box.v4.x, box.v4.y, box.v4.z, 1.0 };
	float ver5[4] = { box.v5.x, box.v5.y, box.v5.z, 1.0 };
	float ver6[4] = { box.v6.x, box.v6.y, box.v6.z, 1.0 };
	float ver7[4] = { box.v7.x, box.v7.y, box.v7.z, 1.0 };
	float ver8[4] = { box.v8.x, box.v8.y, box.v8.z, 1.0 };

	float y_value[8];
	float x_value[8];
	float z_value[8];
	float screen1[3] = { 0.0,0.0,0.0 };
	float screen2[3] = { 0.0,0.0,0.0 };
	float screen3[3] = { 0.0,0.0,0.0 };
	float screen4[3] = { 0.0,0.0,0.0 };
	float screen5[3] = { 0.0,0.0,0.0 };
	float screen6[3] = { 0.0,0.0,0.0 };
	float screen7[3] = { 0.0,0.0,0.0 };
	float screen8[3] = { 0.0,0.0,0.0 };
	transform(mat, ver1, screen1);
	transform(mat, ver2, screen2);
	transform(mat, ver3, screen3);
	transform(mat, ver4, screen4);
	transform(mat, ver5, screen5);
	transform(mat, ver6, screen6);
	transform(mat, ver7, screen7);
	transform(mat, ver8, screen8);
	//X value
	x_value[0] = screen1[0];
	x_value[1] = screen2[0];
	x_value[2] = screen3[0];
	x_value[3] = screen4[0];
	x_value[4] = screen5[0];
	x_value[5] = screen6[0];
	x_value[6] = screen7[0];
	x_value[7] = screen8[0];
	//Y value
	y_value[0] = screen1[1];
	y_value[1] = screen2[1];
	y_value[2] = screen3[1];
	y_value[3] = screen4[1];
	y_value[4] = screen5[1];
	y_value[5] = screen6[1];
	y_value[6] = screen7[1];
	y_value[7] = screen8[1];
	//Z value
	z_value[0] = screen1[2];
	z_value[1] = screen2[2];
	z_value[2] = screen3[2];
	z_value[3] = screen4[2];
	z_value[4] = screen5[2];
	z_value[5] = screen6[2];
	z_value[6] = screen7[2];
	z_value[7] = screen8[2];

	float z_max = z_value[0];
	float z_min = z_value[0];
	float x_max = x_value[0];
	float x_min = x_value[0];
	float y_max = y_value[0];
	float y_min = y_value[0];

	for (int i = 1; i < 8; i++) {
		z_max = max(z_max, z_value[i]);
		z_min = min(z_min, z_value[i]);
		x_max = max(x_max, x_value[i]);
		x_min = min(x_min, x_value[i]);
		y_max = max(y_max, y_value[i]);
		y_min = min(y_min, y_value[i]);
	}

	float result[6] = { x_min, x_max, y_min, y_max, z_min, z_max };
	return result;
}


void Cloud::isCloudIntervel(GzMatrix matrix, float* interval)
{
	boundary box;

	if (object_type == SPHERE)
	{

		box.v1.x = center_x - sphereCloud.radius;
		box.v1.y = center_y + sphereCloud.radius;
		box.v1.z = center_z - sphereCloud.radius;

		box.v2.x = center_x + sphereCloud.radius;
		box.v2.y = center_y + sphereCloud.radius;
		box.v2.z = center_z - sphereCloud.radius;

		box.v3.x = center_x - sphereCloud.radius;
		box.v3.y = center_y - sphereCloud.radius;
		box.v3.z = center_z - sphereCloud.radius;

		box.v4.x = center_x + sphereCloud.radius;
		box.v4.y = center_y - sphereCloud.radius;
		box.v4.z = center_z - sphereCloud.radius;

		box.v5.x = center_x - sphereCloud.radius;
		box.v5.y = center_y + sphereCloud.radius;
		box.v5.z = center_z + sphereCloud.radius;

		box.v6.x = center_x + sphereCloud.radius;
		box.v6.y = center_y + sphereCloud.radius;
		box.v6.z = center_z + sphereCloud.radius;

		box.v7.x = center_x - sphereCloud.radius;
		box.v7.y = center_y - sphereCloud.radius;
		box.v7.z = center_z + sphereCloud.radius;

		box.v8.x = center_x + sphereCloud.radius;
		box.v8.y = center_y - sphereCloud.radius;
		box.v8.z = center_z + sphereCloud.radius;

		float* res = calculateZInterval(matrix, box);
		interval[0] = res[0];
		interval[1] = res[1];
		interval[2] = res[2];
		interval[3] = res[3];
		interval[4] = res[4];
		interval[5] = res[5];

	}

	else if (object_type == CYLINDER) {


		float half_length = cylinderCloud.height / 2.0;

		box.v1.x = center_x - half_length;
		box.v1.y = center_y + cylinderCloud.radius;
		box.v1.z = center_z - cylinderCloud.radius;

		box.v2.x = center_x + half_length;
		box.v2.y = center_y + cylinderCloud.radius;
		box.v2.z = center_z - cylinderCloud.radius;

		box.v3.x = center_x - half_length;
		box.v3.y = center_y - cylinderCloud.radius;
		box.v3.z = center_z - cylinderCloud.radius;

		box.v4.x = center_x + half_length;
		box.v4.y = center_y - cylinderCloud.radius;
		box.v4.z = center_z - cylinderCloud.radius;

		box.v5.x = center_x - half_length;
		box.v5.y = center_y + cylinderCloud.radius;
		box.v5.z = center_z + cylinderCloud.radius;

		box.v6.x = center_x + half_length;
		box.v6.y = center_y + cylinderCloud.radius;
		box.v6.z = center_z + cylinderCloud.radius;

		box.v7.x = center_x - half_length;
		box.v7.y = center_y - cylinderCloud.radius;
		box.v7.z = center_z + cylinderCloud.radius;

		box.v8.x = center_x + half_length;
		box.v8.y = center_y - cylinderCloud.radius;
		box.v8.z = center_z + cylinderCloud.radius;

		float* res = calculateZInterval(matrix, box);
		interval[0] = res[0];
		interval[1] = res[1];
		interval[2] = res[2];
		interval[3] = res[3];
		interval[4] = res[4];
		interval[5] = res[5];
	}

	else if (object_type == 3) {
		float half_length = blockCloud.length / 2.0;
		float half_width = blockCloud.width / 2.0;
		float half_height = blockCloud.height / 2.0;

		box.v1.x = center_x - half_length;
		box.v1.y = center_y + half_width;
		box.v1.z = center_z - half_height;

		box.v2.x = center_x + half_length;
		box.v2.y = center_y + half_width;
		box.v2.z = center_z - half_height;

		box.v3.x = center_x - half_length;
		box.v3.y = center_y - half_width;
		box.v3.z = center_z - half_height;

		box.v4.x = center_x + half_length;
		box.v4.y = center_y - half_width;
		box.v4.z = center_z - half_height;;

		box.v5.x = center_x - half_length;
		box.v5.y = center_y + half_width;
		box.v5.z = center_z + half_height;

		box.v6.x = center_x + half_length;
		box.v6.y = center_y + half_width;
		box.v6.z = center_z + half_height;

		box.v7.x = center_x - half_length;
		box.v7.y = center_y - half_width;
		box.v7.z = center_z + half_height;

		box.v8.x = center_x + half_length;
		box.v8.y = center_y - half_width;
		box.v8.z = center_z + half_height;

		float* res = calculateZInterval(matrix, box);
		interval[0] = res[0];
		interval[1] = res[1];
		interval[2] = res[2];
		interval[3] = res[3];
		interval[4] = res[4];
		interval[5] = res[5];
	}

	else if (object_type == 4)
	{
		float half_length = cubeCloud.length / 2.0;

		box.v1.x = center_x - half_length;
		box.v1.y = center_y + half_length;
		box.v1.z = center_z - half_length;

		box.v2.x = center_x + half_length;
		box.v2.y = center_y + half_length;
		box.v2.z = center_z - half_length;

		box.v3.x = center_x - half_length;
		box.v3.y = center_y - half_length;
		box.v3.z = center_z - half_length;

		box.v4.x = center_x + half_length;
		box.v4.y = center_y - half_length;
		box.v4.z = center_z - half_length;

		box.v5.x = center_x - half_length;
		box.v5.y = center_y + half_length;
		box.v5.z = center_z + half_length;

		box.v6.x = center_x + half_length;
		box.v6.y = center_y + half_length;
		box.v6.z = center_z + half_length;

		box.v7.x = center_x - half_length;
		box.v7.y = center_y - half_length;
		box.v7.z = center_z + half_length;

		box.v8.x = center_x + half_length;
		box.v8.y = center_y - half_length;
		box.v8.z = center_z + half_length;

		float* res = calculateZInterval(matrix, box);
		interval[0] = res[0];
		interval[1] = res[1];
		interval[2] = res[2];
		interval[3] = res[3];
		interval[4] = res[4];
		interval[5] = res[5];
	}
}