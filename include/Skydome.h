#pragma once
#include <unordered_map>

class Skydome
{
public:
	Skydome();
	~Skydome();

	float3 GetSkyColor(Ray& ray);

private:

	int width, height;
	float* pixels;
};