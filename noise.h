#pragma once
#include "Gz.h";
#include "object.h"
#define MOD3 vec3(.1031,.11369,.13787)
class Noise { // Define a noise api

public:

	//Variable declartion 
	const float F3 = 0.3333333;
	const float G3 = 0.1666667;
	float rot1[3][3] = { -0.37, 0.36, 0.85,-0.14,-0.93, 0.34,0.92, 0.01,0.4 };
	float rot2[3][3] = { -0.55,-0.39, 0.74, 0.33,-0.91,-0.24,0.77, 0.12,0.63 };
	float rot3[3][3] = { -0.71, 0.52,-0.47,-0.08,-0.72,-0.68,-0.7,-0.45,0.56 };

	//Constructor
	Noise();
	~Noise();


	Vec3f random3(Vec3f c);
	float simplex3d(Vec3f p);
	float simplex3d_fractal(Vec3f m);
};