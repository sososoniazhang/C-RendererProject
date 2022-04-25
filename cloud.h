
#pragma once
#include	"gz.h"
#include  "object.h"
#include "noise.h"
/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100	/* how many matrix pushes allowed */
#define SPHERE 1
#define CYLINDER 2
#define BLOCK 3
#define CUBE 4

struct boundary
{
	vertex v1;
	vertex v2;
	vertex v3;
	vertex v4;
	vertex v5;
	vertex v6;
	vertex v7;
	vertex v8;
};

class Cloud {			/* define a Cloud Api*/


public:
	//Variable declaration

	//Cloud Center
	float center_x;
	float center_y;
	float center_z;

	sphere sphereCloud;
	cylinder cylinderCloud;
	block blockCloud;
	cube cubeCloud;

	int object_type = 0;
	//Transform matrix
	GzMatrix Xcloud[MATLEVELS];

	short	matlevel;	/* top of stack - current xform */
	GzCamera c_camera;

	Noise* noise = new Noise();

	// Constructors
	Cloud(float centerX, float centerY, float centerZ);
	~Cloud();

	//Default Function declaration
	//int CloudDefault();

	//Update Function declaration
	int putSphere(float radius);
	int putCylinder(float radius, float height);
	int putBlock(float length, float width, float height);
	int putCube(float length);

	// Function declaration
	float isCloud(GzCoord point); //Check if there is cloud on this point
	float isCloudPerlin(float point[4]);
	void isCloudIntervel(GzMatrix matrix, float* interval);// This function will find where the light get in and where it get out, it will save the z min and z max in interval pointer
};