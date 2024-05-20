#pragma once

class Plane
{
public:
	Plane(const float3& origin, const float3& normal) : origin(origin), normal(normal) {};

	bool Intersect(Ray& ray) const;

	inline float3 GetNormal() const
	{
		return normal;
	}

	float3 normal = float3(0, 1, 0);

	float3 origin = float3(0);

	uint8_t material = 0;

	float widht = 2.0f;
	
	float height = 2.0f;
};