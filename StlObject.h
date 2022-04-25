#pragma once
#include "object.h"

#ifndef GzTexture
#define GzTexture	GzPointer
#endif

//struct vertex {
//	float x;
//	float y;
//	float z;
//};

class StlObject {			/* define a renderer */
public:
	char* fname;
	char* buffer;
	float MAX_X;
	float MIN_X;
	float MAX_Y;
	float MIN_Y;
	float MAX_Z;
	float MIN_Z;
	char* iniptr;
	int nTri;

	// Constructors
	StlObject(char* filename);
	~StlObject();

	int readSTL();
	int getTri(int triID, vertex* v1, vertex* v2, vertex* v3, vertex* norm);
	int n_getTri(int triID, vertex* v1, vertex* v2, vertex* v3, vertex norm1, vertex norm2, vertex norm3);

private:
	vertex norm_intep(vertex point);
};