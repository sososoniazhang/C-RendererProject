#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"StlObject.h"
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;

vertex operator ^ (const vertex&, const vertex&);
vertex operator - (const vertex&, const vertex&);
bool operator == (const vertex&, const vertex&);
bool pt_match(vertex input, vertex a, vertex b, vertex c);
void unit(vertex&);
bool isBinarySTL(char*);



StlObject::StlObject(char* filename) {
	MAX_X = 0;
	MIN_X = 0;
	MAX_Y = 0;
	MIN_Y = 0;
	MAX_Z = 0;
	MIN_Z = 0;
	buffer = NULL;
	iniptr = NULL;
	nTri = 0;
	fname = (char*) malloc((strlen(filename)+2));
	if(fname == NULL) return;
	strcpy(fname, filename);
	//memcpy((void*)fname, (void*)filename,28);
}

StlObject::~StlObject() {
	free((void*)fname);
	if (buffer != NULL) {
		free((void*)buffer);
	}
}

int StlObject::readSTL() {
	ifstream ifs(fname, ifstream::binary);
	filebuf* pbuf = ifs.rdbuf();

	auto size = pbuf->pubseekoff(0, ifs.end);
	pbuf->pubseekpos(0);
	buffer = new char[(size_t)size];
	pbuf->sgetn(buffer, size);
	if (!isBinarySTL(buffer)) return -1;

	iniptr = buffer;
	iniptr += 80;  // Skip past the header.
	nTri = *(int*)(iniptr);
	iniptr += 4;   // Skip past the number of triangles.
	return 0;
}

int StlObject::getTri(int triID, vertex* v1, vertex* v2, vertex* v3, vertex* norm) {
	char* triptr = iniptr + (triID * (12 * 4 + 2));
	norm->x = *(float*)(triptr);
	norm->y = *(float*)(triptr + 4);
	norm->z = *(float*)(triptr + 8);
	triptr += 12;

	v1->x = *(float*)(triptr);
	v1->y = *(float*)(triptr + 4);
	v1->z = *(float*)(triptr + 8);
	triptr += 12;

	v2->x = *(float*)(triptr);
	v2->y = *(float*)(triptr + 4);
	v2->z = *(float*)(triptr + 8);
	triptr += 12;

	v3->x = *(float*)(triptr);
	v3->y = *(float*)(triptr + 4);
	v3->z = *(float*)(triptr + 8);
	triptr += 12;

	const float eps = (float)1.0e-9;

	if (abs(norm->x) < eps && abs(norm->y) < eps && abs(norm->z) < eps) {
		vertex u, v;
		u = *v2 - *v1;
		v = *v3 - *v1;
		*norm = u ^ v;
		unit(*norm);
	}

	return 0;
}

int StlObject::n_getTri(int triID, vertex* v1, vertex* v2, vertex* v3, vertex norm1, vertex norm2, vertex norm3) {
	char* triptr = iniptr + (triID * (12 * 4 + 2));
	triptr += 12;

	v1->x = *(float*)(triptr);
	v1->y = *(float*)(triptr + 4);
	v1->z = *(float*)(triptr + 8);
	triptr += 12;

	v2->x = *(float*)(triptr);
	v2->y = *(float*)(triptr + 4);
	v2->z = *(float*)(triptr + 8);
	triptr += 12;

	v3->x = *(float*)(triptr);
	v3->y = *(float*)(triptr + 4);
	v3->z = *(float*)(triptr + 8);
	triptr += 12;

	norm1 = norm_intep(*v1);
	norm2 = norm_intep(*v2);
	norm3 = norm_intep(*v3);

	return 0;
}

vertex StlObject::norm_intep(vertex point) {
	vertex norm, result;
	vertex v1, v2, v3;
	result.x = 0;
	result.y = 0;
	result.z = 0;
	const float eps = (float)1.0e-9;

	for (int i = 0; i < nTri; i++) {
		char* triptr = iniptr + (i * (12 * 4 + 2));
		norm.x = *(float*)(triptr);
		norm.y = *(float*)(triptr + 4);
		norm.z = *(float*)(triptr + 8);
		triptr += 12;

		v1.x = *(float*)(triptr);
		v1.y = *(float*)(triptr + 4);
		v1.z = *(float*)(triptr + 8);
		triptr += 12;

		v2.x = *(float*)(triptr);
		v2.y = *(float*)(triptr + 4);
		v2.z = *(float*)(triptr + 8);
		triptr += 12;

		v3.x = *(float*)(triptr);
		v3.y = *(float*)(triptr + 4);
		v3.z = *(float*)(triptr + 8);
		triptr += 12;

		if (abs(norm.x) < eps && abs(norm.y) < eps && abs(norm.z) < eps) {
			vertex u, v;
			u = v2 - v1;
			v = v3 - v1;
			norm = u ^ v;
			unit(norm);
		}

		if(pt_match(point, v1, v2, v3)){
			result.x += norm.x;
			result.y += norm.y;
			result.z += norm.z;
		}
	}
	unit(result);
	return result;
}


bool pt_match(vertex input, vertex a, vertex b, vertex c) {
	return ((input == a) || (input == b) || (input == c));
}

vertex operator ^ (const vertex& a, const vertex& b) {
	// Cross product.
	vertex result;
	result.x = a.y * b.z - a.z * b.y;
	result.y = a.z * b.x - a.x * b.z;
	result.z = a.x * b.y - a.y * b.x;
	return(result);
}

vertex operator - (const vertex& a, const vertex& b) {
	// Subtraction.
	vertex result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	return(result);
}

bool operator == (const vertex& a, const vertex& b) {
	// equal.
	const float eps = 1;
	return ((abs(a.x - b.x) < eps) && (abs(a.y - b.y) < eps) && (abs(a.z - b.z) < eps));
}

void unit(vertex& v) {
	// Normalize a vector.
	float vmod = pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2);
	vmod = sqrt(vmod);

	if (vmod > (float)1.0e-9) {
		v.x /= vmod;
		v.y /= vmod;
		v.z /= vmod;
	}
}



bool isBinarySTL(char* buffer) {

	// Determine if a file is binary STL.

	bool bbinary = true;
	size_t spnsz, spnsz0;

	// Look for the first non-space character.

	spnsz = strspn(buffer, " ");

	char ctstr[6];  // Enough space to hold "solid\0" and "facet\0".

	// Copy the first five characters from the location of the first non-space
	// character to ctstr.

	strncpy_s(ctstr, &buffer[spnsz], 5);

	ctstr[5] = '\0';
	char csolid[] = "solid\0";

	// If this is an ASCII STL file, then the first string should be "solid".

	if (!strcmp(ctstr, csolid)) {
		// This file might be binary or text. To be certain we need to do a further test.
		// Read past the next new line. If this is a text file, there should be a newline.
		// The next token should be 'facet'.

		spnsz0 = 5 + spnsz;

		char* pch = strchr(&buffer[spnsz0], '\n');  // Look for the first instance of '\n'.
		// If pch is NULL then this is a binary STL file.
		if (pch) {
			pch++;

			spnsz = strspn(pch, " "); // Look for the first instance not of ' '.
			spnsz0 = spnsz;

			spnsz = strcspn(pch + spnsz0, " "); // Look for the first instance of ' '.

			if (spnsz == 5) {
				// Check for 'facet'.
				strncpy_s(ctstr, pch + spnsz0, 5);
				ctstr[5] = '\0';

				char cfacet[] = "facet\0";
				if (!strcmp(ctstr, cfacet)) {
					// This file is beyond reasonable doubt ASCII STL.
					bbinary = false;
				}
			}
		}
	}

	return(bbinary);
}
