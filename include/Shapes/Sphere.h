#pragma once

class Sphere
{
public:
	Sphere() = default;
	Sphere(const float3& center, float radius) : center(center), radius(radius) {};
	void SphereExit(Ray& ray) const;

	float3 center = float3(0, 0, 0);
	float radius = 0.5f;
	uint8_t material = 0;
};