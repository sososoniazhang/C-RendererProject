#pragma once
#include "StlObject.h"
#include <limits>

/* defaults */
#define	MAXLAYER 100
/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

class GzDepthMap {			/* define a depth map */
public:
	unsigned int	xres; //this is the x range of the map
	unsigned int	yres; // this is the y range of the map

	float* zbuffer; // array saves the shadow z value
	// ---- for 3d ------------------------------
	bool threeD;
	float* zbuffers[MAXLAYER];
	int layer; // this is for 3d z buffer index 0 --- (layer-1)
	float zlow;
	float zhigh;
	int step;
	// zhigh = zlow + (layer-1) * step
	// ---- for 3d -----------------------------
	GzCamera lightSource;
	GzColor lightColor;

	short		    matlevel;	        /* top of stack - current xform */
	GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
	GzMatrix		Xsp;
	StlObject* objptr; // this stores the obj pointers to generate the dpeth map
	//float depth; // ************ this parameter is set for the ray marching, for transform matrix NOT NEEDED NOW

	float maxZ; // this is the max z valur in light space

	// for rendering the map
	GzPixel* pixelbuffer;		/* frame buffer array */
	char* framebuffer;




	// Constructors
	GzDepthMap(int zmax = INT_MAX, bool threeD = false);
	~GzDepthMap();

	int startUp();

	// put light source
	int GzPutLightSource(GzCamera lightSource);//done
	int GzPushMatrix(GzMatrix matrix); //done
	int GzPopMatrix(); //done

	// storing the vertex coord
	float v1z, v2z, v3z, v1x, v2x, v3x, v1y, v2y, v3y;
	/***********************************************/
/* distance from view plane to focal point */
	float camDis;
	// for storing A, B, C, D for z interpolation
	float A, B, C, D;



	// generate depth map
	int GzCreateDepthMapFromSTL(StlObject* objptr);
	int GzCreateDepthMapFromTri(int	numParts, GzToken* nameList, GzPointer* valueList);
	float GzFindZ(GzCoord position);

	int GzFlushDisplay2File(FILE* outfile);
	int GzFlushDisplay2FrameBuffer();
	void getShadowMap();
	int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
	int GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z);
	void putColor(GzColor color);

	// Object Translation
	int GzRotXMat(float degree, GzMatrix mat);
	int GzRotYMat(float degree, GzMatrix mat);
	int GzRotZMat(float degree, GzMatrix mat);
	int GzTrxMat(GzCoord translate, GzMatrix mat);
	int GzScaleMat(GzCoord scale, GzMatrix mat);



	// Extra methods: NOT part of API - just for general assistance */
	inline int ARRAY(int x, int y) { return (x + y * xres); }	/* simplify fbuf indexing */
	inline short	ctoi(float color) { return(short)((int)(color * ((1 << 12) - 1))); }		/* convert float color to GzIntensity short */
	float determinate(float a, float b, float c, float d) { return a * d - c * b; }
	float interpolateZ(int x, int y, float* a, float* b, float* c, float* d) { return -(*a * x + *b * y + *d) / *c; }
	float degToRadius(float degree) { return degree * 3.1415926 / 180; }
	float dot(GzCoord x, GzCoord y) { return (*x) * (*y) + (*(x + 1)) * (*(y + 1)) + (*(x + 2)) * (*(y + 2)); }
	bool inBound(int i, int j, int xres, int yres) {
		if (i < 0 || i >= xres || j < 0 || j >= yres) {
			return false;
		}
		return true;
	}

	//helper functions:
	int transWS(vertex* position);
	int readTriangle(vertex a, vertex b, vertex c);
	int drawMap();
	float interpolateDepth2d(float x, float y, int layer = 0);
	float interpolateDepth3d(float x, float y, float z);
	void swapCord(GzCoord* vertexList, int a, int b);
	void getCoordTri(GzCoord* vertexList);
	int sortVert(GzCoord* vertexList, GzCoord* normalList);
	void drawLineCW(GzCoord* vertexList);
	int pointInLine(float lineParams[], int x, int y);
	void crossProd(float* v1, float* v2, float* a, float* b, float* c);
	void findABCD(float* a, float* b, float* c, float* d);
	void fourDto3d(float* newV);
	void project(GzMatrix transX, GzCoord oldV, GzCoord newV, int num);
	float generalLinearInt(float z1, float z2, float z3, int x, int y);
	void normalize(GzCoord x);
	void getXiwZ(GzCoord* z, GzCoord c, GzCoord l);
	void getXiwY(GzCoord* y, GzCoord up, GzCoord z);
	void getXiwX(GzCoord* x, GzCoord y, GzCoord z);
	void getXiw(GzCamera* cameraPtr);
};