#pragma once

struct Mesh
{
	Mesh(const char* filename);
	Triangle* triangles; // triangle array
	TriangleNormal* normals; // normal array
	int numTriangles = 0; // number of triangles
};