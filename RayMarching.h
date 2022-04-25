#pragma once

//#include "stdafx.h"
//#include "CS580HW.h"
//#include "ApplicationCloud.h"
#include "Gz.h"
#include "cloud.h"
#include "object.h"

#define XMIN 0 
#define XMAX 1 
#define YMIN 2 
#define YMAX 3 
#define ZMIN 4 
#define ZMAX 5 


struct ShadeBuffer {
	float* Solidbuffer;
	float* CloudRangebuffer;
	float* CloudCoefbuffer[100];
};

class RayMarch {
	Cloud* cloud;
	GzMatrix Xsw;

public:
	unsigned short	xres;
	unsigned short	yres;
	GzPixel* pixelbuffer;		/* frame buffer array */
	char* framebuffer;
	GzCamera		rm_camera;
	GzMatrix	inv_Xiw;
	GzMatrix	Xsp;
	float boundary[6];
	float zboundary[6];

	short		    matlevel;	        /* top of stack - current xform */
	GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */

	// Constructors
	RayMarch(int xRes, int yRes, GzCamera camera, Cloud* new_cloud, GzPixel* pbuffer, char* fbuffer);
	~RayMarch();

	//Methods
	int camera_change(GzCamera camera);
	float RayMarching(GzCoord pt_s);
	
	void RMrender();
	void RMrender_buffer(ShadeBuffer* Sbuffer);

	int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
	int FlushDisplay2FrameBuffer();
	int RmDefault(int mred, int mgreen, int mblue);

	int PushMatrix(GzMatrix	matrix);
	int PopMatrix();

	inline int ARRAY(int x, int y) { return (x + y * xres); }	/* simplify fbuf indexing */
};