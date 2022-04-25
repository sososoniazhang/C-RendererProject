// Application4.h: interface for the Application4 class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APPLICATIONC_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
#define AFX_APPLICATION3_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Application.h"
#include "StlObject.h"
#include "cloud.h"
#include "RayMarching.h"

class ApplicationCloud : public Application
{
	StlObject* model;
	GzPixel* app_pixelbuffer;		/* frame buffer array */
	char* app_framebuffer;
	GzDepthMap* lightShader;

public:
	RayMarch* rmRender[10];
	Cloud* cloud_obj[10];

	ApplicationCloud();
	virtual ~ApplicationCloud();

	int	Initialize();
	virtual int Render();
	int Clean();
};

#endif // !defined(AFX_APPLICATION4_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
