#include "stdafx.h"
#include "RayMarching.h"
//#include "Matrix.h"

#define CLOUD_REF 0.005

#define Freqency 200

#define PI (float) 3.14159265358979323846

extern bool inverseMatrix(GzMatrix mat, GzMatrix inverse);

//////////////////////////////////////Helpers//////////////////////////////////////
static void xform(GzCoord input, GzMatrix xform_matrix, GzCoord result) {
	float temp[4] = { 0,0,0,0 };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (j == 3) {
				temp[i] = temp[i] + (xform_matrix[i][j] * 1);
			}
			else {
				temp[i] = temp[i] + (xform_matrix[i][j] * input[j]);
			}

		}
	}
	for (int i = 0; i < 3; i++) {
		result[i] = temp[i] / temp[3];
	}
}

static void crossP(GzMatrix matrix_a, GzMatrix matrix_b, GzMatrix result) {
	memset((void*)result, 0, sizeof(GzMatrix));
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				result[i][j] = result[i][j] + (matrix_a[i][k] * matrix_b[k][j]);
			}
		}
	}
}

//////////////////////////////////////Main//////////////////////////////////////
RayMarch::RayMarch(int xRes, int yRes, GzCamera camera, Cloud* new_cloud, GzPixel* pbuffer, char* fbuffer) {
	xres = xRes;
	yres = yRes; 

	//framebuffer = (char*)malloc((size_t)3 * xres * yres);
	//pixelbuffer = (GzPixel*)malloc((size_t)(sizeof(GzPixel) * xres * yres));
	pixelbuffer = pbuffer;
	framebuffer = fbuffer;
	matlevel = -1;

	memset((void*)Xsw, 0, sizeof(GzMatrix));

	camera_change(camera);
	cloud = new_cloud;
	
}

RayMarch::~RayMarch() {
	//free((void*)framebuffer);
	//free((void*)pixelbuffer);
}

//Methods
int RayMarch::camera_change(GzCamera camera) {
	rm_camera = camera;

	memset((void*)inv_Xiw, 0, sizeof(GzMatrix));
	/*************  Setting Xiw  **************/
	GzCoord c_l, up_p, x_vec, y_vec, z_vec;
	for (int i = 0; i < 3; i++) {
		c_l[i] = (rm_camera.lookat[i] - rm_camera.position[i]);
	}
	float c_l_norm = sqrt((c_l[0] * c_l[0]) + (c_l[1] * c_l[1]) + (c_l[2] * c_l[2]));
	for (int i = 0; i < 3; i++) {
		z_vec[i] = c_l[i] / c_l_norm;
	}

	float up_z = (rm_camera.worldup[0] * z_vec[0]) + (rm_camera.worldup[1] * z_vec[1]) + (rm_camera.worldup[2] * z_vec[2]);
	for (int i = 0; i < 3; i++) {
		up_p[i] = rm_camera.worldup[i] - up_z * z_vec[i];
	}
	float up_p_norm = sqrt((up_p[0] * up_p[0]) + (up_p[1] * up_p[1]) + (up_p[2] * up_p[2]));
	for (int i = 0; i < 3; i++) {
		y_vec[i] = up_p[i] / up_p_norm;
	}

	x_vec[0] = (y_vec[1] * z_vec[2]) - (y_vec[2] * z_vec[1]);
	x_vec[1] = -((y_vec[0] * z_vec[2]) - (y_vec[2] * z_vec[0]));
	x_vec[2] = (y_vec[0] * z_vec[1]) - (y_vec[1] * z_vec[0]);

	memset((void*)&rm_camera.Xiw, 0, sizeof(GzMatrix));
	rm_camera.Xiw[0][0] = x_vec[0];
	rm_camera.Xiw[0][1] = x_vec[1];
	rm_camera.Xiw[0][2] = x_vec[2];
	rm_camera.Xiw[1][0] = y_vec[0];
	rm_camera.Xiw[1][1] = y_vec[1];
	rm_camera.Xiw[1][2] = y_vec[2];
	rm_camera.Xiw[2][0] = z_vec[0];
	rm_camera.Xiw[2][1] = z_vec[1];
	rm_camera.Xiw[2][2] = z_vec[2];

	rm_camera.Xiw[0][3] = -((x_vec[0] * rm_camera.position[0]) + (x_vec[1] * rm_camera.position[1]) + (x_vec[2] * rm_camera.position[2]));
	rm_camera.Xiw[1][3] = -((y_vec[0] * rm_camera.position[0]) + (y_vec[1] * rm_camera.position[1]) + (y_vec[2] * rm_camera.position[2]));
	rm_camera.Xiw[2][3] = -((z_vec[0] * rm_camera.position[0]) + (z_vec[1] * rm_camera.position[1]) + (z_vec[2] * rm_camera.position[2]));

	rm_camera.Xiw[3][3] = 1;

	if (!inverseMatrix(rm_camera.Xiw, inv_Xiw)) {
		exit(-10);
	}

	/*************  Setting Xpi  **************/
	memset((void*)&rm_camera.Xpi, 0, sizeof(GzMatrix));
	float inv_d = tan(rm_camera.FOV * PI / 360);
	rm_camera.Xpi[0][0] = 1;
	rm_camera.Xpi[1][1] = 1;
	rm_camera.Xpi[2][2] = inv_d;
	rm_camera.Xpi[3][3] = 1;
	rm_camera.Xpi[3][2] = inv_d;

	/*************  Setting Xsp  **************/
	memset((void*)Xsp, 0, sizeof(GzMatrix));
	Xsp[0][0] = yres / 2;
	Xsp[1][1] = -yres / 2;
	Xsp[2][2] = MAXINT;
	Xsp[3][3] = 1;
	Xsp[0][3] = xres / 2;
	Xsp[1][3] = yres / 2;

	PushMatrix(Xsp);
	PushMatrix(rm_camera.Xpi);
	PushMatrix(rm_camera.Xiw);
	return 0;
}

float RayMarch::RayMarching(GzCoord pt_s) {
	float step = (zboundary[ZMAX] - zboundary[ZMIN]) / Freqency;
	float inv_d = tan(rm_camera.FOV * PI / 360);
	float cloud_coef = 0.0;
	GzCoord march_point,world_point;
	float result = 0;
	for (int i = 0; i < 3; i++) {
		march_point[i] = pt_s[i];
	}
	for (int i = 0; i < Freqency; i++) {
		march_point[Z] = zboundary[ZMIN] + (i * step);
		march_point[X] = pt_s[X] * ((march_point[Z] * inv_d) + 1);
		march_point[Y] = pt_s[Y] * ((march_point[Z] * inv_d) + 1);
		xform(march_point, inv_Xiw, world_point);
		//cloud_coef = (rand() % 2) * 1; cloud_coef *
		result += cloud->isCloud(world_point) * CLOUD_REF /*(1 / Freqency) * 100*/;
	}
	return result;
}

void RayMarch::RMrender() {
	GzCoord pixel_point;
	float color = 0;
	int c12b = 0;
	float xstep = 2.0 / yres;
	float ystep = 2.0 / yres;
	cloud->isCloudIntervel(Ximage[matlevel], boundary);
	cloud->isCloudIntervel(rm_camera.Xiw, zboundary);
	int xstart = boundary[XMIN] - 1;
	int xend = boundary[XMAX] + 1;
	int ystart = boundary[YMIN] - 1;
	int yend = boundary[YMAX] + 1;
	if (xstart < 0) xstart = 0;
	if (xend > xres - 1) xend = xres - 1;
	if (xstart > xend) return;
	if (ystart < 0) ystart = 0;
	if (yend > yres - 1) yend = yres - 1;
	if (ystart > yend) return;
	for (int i = xstart; i <= xend; i++) {
		for (int j = ystart; j <= yend; j++) {
			float xcoord = -(xres / 2 * xstep) + (((double)i) * xstep);
			float ycoord = -(yres / 2 * ystep) + (((double)j) * ystep);
			pixel_point[X] = xcoord;
			pixel_point[Y] = -ycoord;
			pixel_point[Z] = 0;
			color = RayMarching(pixel_point);
			if (color > 1) c12b = 4095;
			else if (color < 0) c12b = 0;
			else c12b = (int)(4095.0 * color);
			GzPut(i, j, 4095, 4095, 4095, c12b, zboundary[ZMIN]);
		}
	}
}

void RayMarch::RMrender_buffer(ShadeBuffer* Sbuffer) {
	GzCoord pixel_point;
	float color = 0;
	int c12b = 0;
	float xstep = 2.0 / xres;

	float ystep = 2.0 / yres;
	cloud->isCloudIntervel(Ximage[matlevel], boundary);
	cloud->isCloudIntervel(rm_camera.Xiw, zboundary);
	int xstart = boundary[XMIN] - 1;
	int xend = boundary[XMAX] + 1;
	int ystart = boundary[YMIN] - 1;
	int yend = boundary[YMAX] + 1;
	if (xstart < 0) xstart = 0;
	if (xend > xres - 1) xend = xres - 1;
	if (xstart > xend) return;
	if (ystart < 0) ystart = 0;
	if (yend > yres - 1) yend = yres - 1;
	if (ystart > yend) return;
	for (int i = xstart; i <= xend; i++) {
		for (int j = ystart; j <= yend; j++) {
			float xcoord = -1.0 + (((double)i) * xstep);
			float ycoord = -1.0 + (((double)j) * ystep);
			pixel_point[X] = xcoord;
			pixel_point[Y] = -ycoord;
			pixel_point[Z] = 0;
			color = RayMarching(pixel_point);
			if (color > 1) c12b = 4095;
			else if (color < 0) c12b = 0;
			else c12b = (int)(4095.0 * color);
			GzPut(i, j, 4095, 4095, 4095, c12b, zboundary[ZMIN]);
		}
	}
}

int RayMarch::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */
	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		GzPixel* pbuf_ptr = pixelbuffer + ARRAY(i, j);
		if (z < pbuf_ptr->z && a != 0) {
			if (r > 4095) r = 4095;
			if (r < 0) r = 0;
			if (g > 4095) g = 4095;
			if (g < 0) g = 0;
			if (b > 4095) b = 4095;
			if (b < 0) b = 0;
			if (a > 4095) a = 4095;
			if (a < 0) a = 0;
			pbuf_ptr->red = ((pbuf_ptr->red * (4095 - a)) + (r * a)) / 4095;
			pbuf_ptr->green = ((pbuf_ptr->green * (4095 - a)) + (g * a)) / 4095;
			pbuf_ptr->blue = ((pbuf_ptr->blue * (4095 - a)) + (b * a)) / 4095;
			//pbuf_ptr->alpha = a;
			//pbuf_ptr->z = z;
		}
	}
	return GZ_SUCCESS;
}


int RayMarch::FlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	GzPixel* pbuf_ptr = pixelbuffer;
	char* fbuf_ptr = framebuffer;
	for (int j = 0; j < yres; j++) {
		for (int i = 0; i < xres; i++) {
			char wred = (char)((int)(pbuf_ptr->red) >> 4);
			char wgreen = (char)((int)(pbuf_ptr->green) >> 4);
			char wblue = (char)((int)(pbuf_ptr->blue) >> 4);
			fbuf_ptr[3 * ARRAY(i, j)] = wblue;
			fbuf_ptr[3 * ARRAY(i, j) + 1] = wgreen;
			fbuf_ptr[3 * ARRAY(i, j) + 2] = wred;
			pbuf_ptr++;
		}
	}
	return GZ_SUCCESS;
}

int RayMarch::RmDefault(int mred, int mgreen, int mblue)
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */
	GzPixel* pbuf_ptr = pixelbuffer;
	//int mred = 135;
	//int mgreen = 206;
	//int mblue = 235;
	for (int i = 0; i < (xres * yres); i++) {
		pbuf_ptr->red = mred << 4;
		pbuf_ptr->green = mgreen << 4;
		pbuf_ptr->blue = mblue << 4;
		pbuf_ptr->alpha = 4095;
		pbuf_ptr->z = MAXINT;
		pbuf_ptr++;
	}
	return GZ_SUCCESS;
}

int RayMarch::PushMatrix(GzMatrix	matrix)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	GzMatrix rotation_matrix;
	GzMatrix empty_matrix = {
	1.0,	0.0,	0.0,	0.0,
	0.0,	1.0,	0.0,	0.0,
	0.0,	0.0,	1.0,	0.0,
	0.0,	0.0,	0.0,	1.0
	};
	if (matlevel == -1) {
		matlevel = matlevel + 1;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Ximage[matlevel][i][j] = matrix[i][j];
			}
		}
		return GZ_SUCCESS;
	}
	else if (matlevel < 99) {
		matlevel = matlevel + 1;
		crossP(Ximage[matlevel - 1], matrix, Ximage[matlevel]);
		return GZ_SUCCESS;
	}
	else {
		return GZ_FAILURE;
	}
}

int RayMarch::PopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel < 0) {
		return GZ_FAILURE;
	}
	else {
		matlevel = matlevel - 1;
		return GZ_SUCCESS;
	}

}

