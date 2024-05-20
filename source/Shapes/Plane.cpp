#include "precomp.h"
#include "include/Shapes/Plane.h"

bool Plane::Intersect(Ray& ray) const
{
	float3 n = normal;
	float denom = dot(n, ray.D);
	if (denom < 1e-6)
	{
		n = -normal;
	}

	denom = dot(n, ray.D);

	const float3 p0l0 = origin - ray.O;
	const float t = dot(p0l0, n) / denom;
	if (t < 0 || ray.t < t)
	{
		return false;
	}

	/*
	float3 p = ray.O + t * ray.D;

	float3 dir = p - origin;
	float3 normalDir = normalize(dir);

	float3 u = normalize(cross(normalDir, n));
	float3 v = normalize(cross(u, normal));

	float boundu = dot(u, dir);
	float boundv = dot(v, dir);

	float minX = -widht * 0.5f;
	float maxX = widht * 0.5f;
	float minY = -height * 0.5f;
	float maxY = height * 0.5f;

	if (boundu < minX || boundu > maxX || boundv > maxY || boundv < minY)
	{
		return false;
	}
	*/

	if (ray.t > t)
	{
		ray.t = t;
		ray.voxel = material;
		ray.normal = normal;
		ray.type = Primitive::PlanePrimitive;
		return true;
	}
	return false;
}