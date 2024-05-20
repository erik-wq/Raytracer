#include "precomp.h"
#include "include/Shapes/Sphere.h"


void Sphere::SphereExit(Ray& ray) const
{
	const float3 sphereDirection = ray.O - center; // order reversed

	const float a = dot(ray.D, ray.D);
	const float b = 2 * dot(sphereDirection, ray.D);
	const float c = dot(sphereDirection, sphereDirection) - radius * radius;

	const float d = b * b - 4 * a * c;
	if (d < 0.0f)
	{
		return;
	}

	const float discSquared = sqrtf(d);
	const float a2 = 2 * a;
	const float t1 = (-b + discSquared) / a2;
	const float t2 = (-b - discSquared) / a2;

	ray.t = t1;
	if (t1 < 0)
	{
		ray.t = t2;
	}
	ray.normal = center;
	ray.type = Primitive::SpherePrimitive;
}