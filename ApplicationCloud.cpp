// Application4.cpp: implementation of the Application4 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment
 */

#include "stdafx.h"
#include "CS580HW.h"
#include "ApplicationCloud.h"
#include "Gz.h"
#include "rend.h"
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>



#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#define new DEBUG_NEW
#endif


#define OUTFILE "output.ppm"
//#define INFILE  "Bunny-Final.stl"
#define INFILE  "pikachu_1gen_flowalistik.STL"
//#define INFILE  "origamix_crane.stl"
//#define INFILE3  "test.asc"
//#define INFILE  "pot4.asc"
//#define INFILE  "tri.asc"

#define XOFF -6
#define YOFF -4
#define ZOFF 0

void shade(GzCoord norm, GzCoord color);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ApplicationCloud::ApplicationCloud()
{
	model = NULL;
	for (int i = 0; i < 10; i++) {
		rmRender[i] = NULL;
		cloud_obj[i] = NULL;
	}
	
	app_framebuffer = NULL;
	app_pixelbuffer = NULL;
}

ApplicationCloud::~ApplicationCloud()
{
	Clean();
}

int ApplicationCloud::Initialize()
{
	/* to be filled in by the app if it sets camera params */

	GzCamera	camera;
	GzCamera	lightSource;
	int		xRes, yRes;		/* display parameters */

	GzToken		nameListShader[9]; 	/* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int		shaderType, interpStyle;
	float		specpower;
	int		status;

	status = 0;


	model = new StlObject(INFILE);
	model->readSTL();

	/*
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/*
	 * initialize the display and the renderer
	 */

	m_nWidth = 500;	 	// frame buffer and display width
	m_nHeight = 800;	// frame buffer and display height
	//m_nWidth = 500;	 	// frame buffer and display width
	//m_nHeight = 800;	// frame buffer and display height

	m_pRender = new GzRender(m_nWidth, m_nHeight);
	m_pRender->GzDefault();

	

	lightShader = new GzDepthMap();
	m_pFrameBuffer = m_pRender->framebuffer;

	/* Translation matrix */
	GzMatrix	scale =
	{
		0.25,	0.0,	0.0,	0.0,
		0.0,	0.25,	0.0,	0.0,
		0.0,	0.0,	0.25,	0.0,
		0.0,	0.0,	0.0,	1.0
	};

	GzMatrix	rotateX =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};

	GzMatrix	rotateY =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};

#if 1 	/* set up app-defined camera if desired, else use camera defaults */
	camera.position[X] = 15;
	camera.position[Y] = 15;
	camera.position[Z] = 17;

	camera.lookat[X] = 0;
	camera.lookat[Y] = 0;
	camera.lookat[Z] = 7;

	camera.worldup[X] = 0.0;
	camera.worldup[Y] = 0.0;
	camera.worldup[Z] = 1.0;

	camera.FOV = 53.7;              /* degrees */

	status |= m_pRender->GzPutCamera(camera);

	lightSource.position[X] = 0.0;
	lightSource.position[Y] = 20;
	lightSource.position[Z] = 10;

	lightSource.lookat[X] = 0.0;
	lightSource.lookat[Y] = 0.0;
	lightSource.lookat[Z] = 0.0;

	lightSource.worldup[X] = 0.0;
	lightSource.worldup[Y] = 1.0;
	lightSource.worldup[Z] = 0.0;

	lightSource.FOV = 53;

	status |= lightShader->GzPutLightSource(lightSource);
	status |= m_pRender->GzPutShader(lightShader);

	/* Light */
	GzColor	lightColor = { 0.2, 0.5, 0.9 };
	lightShader->putColor(lightColor);
#endif 

	/* Start Renderer */
	lightShader->startUp();
	status |= m_pRender->GzBeginRender();

	/* Light */
	//GzLight	light1 = { {0, 0.7071, -0.7071}, {0.8, 0.8, 0.4} };
	//GzLight	light2 = { {0, -0.7071, -0.7071}, {0.5, 0.5, 0.3} };
	//GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.6, 0.6, 0.2} };
	GzLight	ambientlight = { {0, 0, 0}, {0.8, 0.8, 0.8} };

	/* Material property */
	GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
	GzColor ambientCoefficient = { 0.5, 0.5, 0.0 };
	GzColor diffuseCoefficient = { 0.7, 0.7, 0.7 };

	/*
	  renderer is ready for frame --- define lights and shader at start of frame
	*/

	/*
	 * Tokens associated with light parameters
	 */
	//nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
	//valueListLights[0] = (GzPointer)&light1;
	//nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
	//valueListLights[1] = (GzPointer)&light2;
	//nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
	//valueListLights[2] = (GzPointer)&light3;
	status |= m_pRender->GzPutAttribute(1, nameListLights, valueListLights);

	nameListLights[0] = GZ_AMBIENT_LIGHT;
	valueListLights[0] = (GzPointer)&ambientlight;
	status |= m_pRender->GzPutAttribute(1, nameListLights, valueListLights);

	/*
	 * Tokens associated with shading
	 */
	nameListShader[0] = GZ_DIFFUSE_COEFFICIENT;
	valueListShader[0] = (GzPointer)diffuseCoefficient;

	/*
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode
	*/
	nameListShader[1] = GZ_INTERPOLATE;
#if 0
	interpStyle = GZ_COLOR;         /* Gouraud shading */
#else 
	interpStyle = GZ_NORMAL;       /* Phong shading */
#endif

	valueListShader[1] = (GzPointer)&interpStyle;
	nameListShader[2] = GZ_AMBIENT_COEFFICIENT;
	valueListShader[2] = (GzPointer)ambientCoefficient;
	nameListShader[3] = GZ_SPECULAR_COEFFICIENT;
	valueListShader[3] = (GzPointer)specularCoefficient;
	nameListShader[4] = GZ_DISTRIBUTION_COEFFICIENT;
	specpower = 32;
	valueListShader[4] = (GzPointer)&specpower;

	status |= m_pRender->GzPutAttribute(5, nameListShader, valueListShader);

	status |= m_pRender->GzPushMatrix(scale);
	status |= m_pRender->GzPushMatrix(rotateY);
	status |= m_pRender->GzPushMatrix(rotateX);

	status |= lightShader->GzPushMatrix(scale);
	status |= lightShader->GzPushMatrix(rotateY);
	status |= lightShader->GzPushMatrix(rotateX);

	m_pRender->getXmp();

	if (status) exit(GZ_FAILURE);

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
/*
	GzCamera	camera;

	camera.position[X] = -6;
	camera.position[Y] = -6;
	camera.position[Z] = 6;
	
	camera.lookat[X] = 0;
	camera.lookat[Y] = 0;
	camera.lookat[Z] = 0;
	
	camera.worldup[X] = 0.0;
	camera.worldup[Y] = 0.0;
	camera.worldup[Z] = 1.0;
	
	camera.FOV = 40;   

	GzMatrix	Xsw =
	{
		1.0,	0.0,	0.0,	0.0,
		0.0,	1.0,	0.0,	0.0,
		0.0,	0.0,	1.0,	0.0,
		0.0,	0.0,	0.0,	1.0
	};
	
	m_nWidth = 500;	 	// frame buffer and display width
	m_nHeight = 500;	// frame buffer and display height
	app_pixelbuffer = (GzPixel*)malloc((size_t)(sizeof(GzPixel) * m_nWidth * m_nHeight));
	app_framebuffer = (char*)malloc((size_t)3 * m_nWidth * m_nHeight);

	//for (int i = 0; i < 1; i++) {
	//	cloud_obj[i] = new Cloud(i, i-1, i-1);
	//	if (i % 2 == 0) {
	//		cloud_obj[i]->putSphere(i + 1.0);
	//	}
	//	else {
	//		cloud_obj[i]->putCube(i + 1.0);
	//	}
	//	rmRender[i] = new RayMarch(m_nWidth, m_nHeight, camera, cloud_obj[i], app_pixelbuffer, app_framebuffer);
	//}

	cloud_obj[0] = new Cloud(0.0,0.0,0.0);
	cloud_obj[0]->putBlock(5.0,2.0,2.0);
	rmRender[0] = new RayMarch(m_nWidth, m_nHeight, camera, cloud_obj[0], app_pixelbuffer, app_framebuffer);
	//rmRender[0]->RmDefault(135,206,235);
	rmRender[0]->RmDefault(0, 0, 0);
	m_pFrameBuffer = app_framebuffer;
	return(GZ_SUCCESS);
*/
}

int ApplicationCloud::Render()
{
	/*for (int i = 0; i < 1; i++) {
		rmRender[i]->RMrender();
	}
	rmRender[0]->FlushDisplay2FrameBuffer();
	return(GZ_SUCCESS);*/

	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzToken		nameListColor[3];		/* color type names */
	GzPointer	valueListColor[3];	/* color type rgb pointers */
	GzCoord		vertexList[3];		/* vertex position coordinates */
	GzCoord		normalList[3];		/* vertex normals */
	GzTextureIndex  uvList[3];		/* vertex texture map indices */
	GzColor		color;
	char		dummy[256];
	int		status;


	/* Initialize Display */
	status |= m_pRender->GzDefault();  /* init for new frame */

	/*
	* Tokens associated with triangle vertex values
	*/
	nameListTriangle[0] = GZ_POSITION;
	nameListTriangle[1] = GZ_NORMAL;

	// I/O File open
	/*FILE* infile;
	if( (infile  = fopen( INFILE , "r" )) == NULL )
	{
		 AfxMessageBox(_T("The input file was not opened\n" ));
		 return GZ_FAILURE;
	}*/

	/*FILE* outfile;
	if ((outfile = fopen(OUTFILE, "wb")) == NULL)
	{
		AfxMessageBox(_T("The output file was not opened\n"));
		return GZ_FAILURE;
	}*/
	// Open the file for reading using an input fstream.
	/*int nVertex = 0; // Number of vertices read.
	int nFacet = 0;  // Number of facets read.
	ifstream ifs(INFILE, ifstream::binary);
	filebuf* pbuf = ifs.rdbuf();
	auto size = pbuf->pubseekoff(0, ifs.end);
	pbuf->pubseekpos(0);
	char* buffer = new char[(size_t)size];
	pbuf->sgetn(buffer, size);
	if (!isBinarySTL(buffer)) return -1;
	char* bufptr = buffer;
	bufptr += 80;  // Skip past the header.
	bufptr += 4;   // Skip past the number of triangles.

	*/
	vertex normal;
	vertex n1, n2, n3;
	vertex v1, v2, v3;
	/*
	* Walk through the list of triangles, set color
	* and render each triangle
	*/
	for (int i = 0; i < model->nTri; i++) {
		model->getTri(i, &v1, &v2, &v3, &normal);
		//model->n_getTri(i, &v1, &v2, &v3, n1, n2, n3);

		vertexList[0][0] = v1.x + XOFF;
		vertexList[0][1] = v1.y + YOFF;
		vertexList[0][2] = v1.z + ZOFF;

		vertexList[1][0] = v2.x + XOFF;
		vertexList[1][1] = v2.y + YOFF;
		vertexList[1][2] = v2.z + ZOFF;

		vertexList[2][0] = v3.x + XOFF;
		vertexList[2][1] = v3.y + YOFF;
		vertexList[2][2] = v3.z + ZOFF;

		normalList[0][0] = normal.x;
		normalList[0][1] = normal.y;
		normalList[0][2] = normal.z;

		normalList[1][0] = normal.x;
		normalList[1][1] = normal.y;
		normalList[1][2] = normal.z;

		normalList[2][0] = normal.x;
		normalList[2][1] = normal.y;
		normalList[2][2] = normal.z;

		valueListTriangle[0] = (GzPointer)vertexList;
		valueListTriangle[1] = (GzPointer)normalList;

		lightShader->GzCreateDepthMapFromTri(2, nameListTriangle, valueListTriangle);
	}

	for (int i = 0; i < model->nTri; i++) {
		model->getTri(i, &v1, &v2, &v3, &normal);
		//model->n_getTri(i, &v1, &v2, &v3, n1, n2, n3);

		vertexList[0][0] = v1.x + XOFF;
		vertexList[0][1] = v1.y + YOFF;
		vertexList[0][2] = v1.z + ZOFF;

		vertexList[1][0] = v2.x + XOFF;
		vertexList[1][1] = v2.y + YOFF;
		vertexList[1][2] = v2.z + ZOFF;

		vertexList[2][0] = v3.x + XOFF;
		vertexList[2][1] = v3.y + YOFF;
		vertexList[2][2] = v3.z + ZOFF;

		normalList[0][0] = normal.x;
		normalList[0][1] = normal.y;
		normalList[0][2] = normal.z;

		normalList[1][0] = normal.x;
		normalList[1][1] = normal.y;
		normalList[1][2] = normal.z;

		normalList[2][0] = normal.x;
		normalList[2][1] = normal.y;
		normalList[2][2] = normal.z;
		/*
		* Set up shading attributes for each triangle
		*/
		//shade(normalList[0], color);/* shade based on the norm of vert0 */
		//valueListColor[0] = (GzPointer)color;
		//nameListColor[0] = GZ_RGB_COLOR;
		//m_pRender->GzPutAttribute(1, nameListColor, valueListColor);

		/*
		* Set the value pointers to the first vertex of the
		* triangle, then feed it to the renderer
		*/
		valueListTriangle[0] = (GzPointer)vertexList;
		valueListTriangle[1] = (GzPointer)normalList;

		m_pRender->GzPutTriangle(2, nameListTriangle, valueListTriangle);
	}
	//lightShader->getShadowMap();
	//lightShader->GzFlushDisplay2FrameBuffer();

	//m_pRender->GzFlushDisplay2File(outfile); 	/* write out or update display to file*/
	m_pRender->GzFlushDisplay2FrameBuffer();	// write out or update display to frame buffer

	/*
	 * Close file
	 */

	 //if( fclose( infile ) )
	 //  AfxMessageBox(_T( "The input file was not closed\n" ));

	//if (fclose(outfile))
	//	AfxMessageBox(_T("The output file was not closed\n"));

	if (status)
		return(GZ_FAILURE);
	else
		return(GZ_SUCCESS);
}

int ApplicationCloud::Clean()
{
	/*
	 * Clean up and exit
	 */
	if (model != NULL) {
		delete(model);
	}
	if (m_pRender != NULL) {
		delete(m_pRender);
	}
	if (lightShader != NULL) {
		delete(lightShader);
	}
	/*for (int i = 0; i < 10; i++) {
		if (rmRender[i] != NULL) {
			delete(rmRender[i]);
		}
		if (cloud_obj[i] != NULL) {
			delete(cloud_obj[i]);
		}
	}*/
	if (app_pixelbuffer != NULL) {
		free(app_pixelbuffer);
	}
	if (app_framebuffer != NULL) {
		free(app_framebuffer);
	}
	return(GZ_SUCCESS);
}

/*
This doesn't really belong in the application program, but for this
simplified case of a renderer that doesn't do any shading itself, this
is the easiest place to put it.
*/

void shade(GzCoord norm, GzCoord color)
{
	/*GzCoord	light;
	float		coef;

	light[0] = 0.707f;
	light[1] = 0.5f;
	light[2] = 0.5f;

	coef = light[0]*norm[0] + light[1]*norm[1] + light[2]*norm[2];
	if (coef < 0) 	coef *= -1;

	if (coef > 1.0)	coef = 1.0;
	color[0] = coef*0.95f;
	color[1] = coef*0.65f;
	color[2] = coef*0.88f;*/
	color[0] = 1.0f;
	color[1] = 1.0f;
	color[2] = 1.0f;
}
