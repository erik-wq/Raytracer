#include "precomp.h"
#include "include/BHV.h"

float3 Ray::GetNormal() const
{
	// return normal;
#if TEAPOD
	return normal;
#endif

	// return nomral;

	if (type == Primitive::SpherePrimitive)
	{
		return normalize(IntersectionPoint() - normal);
	}

	else if (type == Primitive::GridPrimitive)
	{
		uint x = (*(uint*)&D.x >> 31);
		uint y = (*(uint*)&D.y >> 31);
		uint z = (*(uint*)&D.z >> 31);

		float3 Dsign = float3((float)x, (float)y, (float)z);

		const float3 I1 = (O + t * D) * WORLDSIZE;
		const float3 fG = fracf(I1);
		const float3 d = min3(fG, 1.0f - fG);
		const float mind = min(min(d.x, d.y), d.z);
		const float3 sign = Dsign * 2 - 1;
		const float3 n = float3(mind == d.x ? sign.x : 0, mind == d.y ? sign.y : 0, mind == d.z ? sign.z : 0);
		return normalize(n);
	}
	else if (type == Primitive::PlanePrimitive)
	{
		return normal;
	}
	return float3(0);
}

Cube::Cube(const float3 pos, const float3 size)
{
	// set cube bounds
	b[0] = pos;
	b[1] = pos + size;
}

float Cube::Intersect(const Ray& ray) const
{
#if 1
	// test if the ray intersects the cube
	const int signx = ray.D.x < 0, signy = ray.D.y < 0, signz = ray.D.z < 0;
	float tmin = (b[signx].x - ray.O.x) * ray.rD.x;
	float tmax = (b[1 - signx].x - ray.O.x) * ray.rD.x;
	const float tymin = (b[signy].y - ray.O.y) * ray.rD.y;
	const float tymax = (b[1 - signy].y - ray.O.y) * ray.rD.y;
	if (tmin > tymax || tymin > tmax) goto miss;
	tmin = max(tmin, tymin), tmax = min(tmax, tymax);
	const float tzmin = (b[signz].z - ray.O.z) * ray.rD.z;
	const float tzmax = (b[1 - signz].z - ray.O.z) * ray.rD.z;
	if (tmin > tzmax || tzmin > tmax) goto miss; // yeah c has 'goto' ;)
	if ((tmin = max(tmin, tzmin)) > 0) return tmin;
miss:
	return 1e34f;
#else
	__m128 sign = _mm_c 

#endif
}

bool Cube::Contains(const float3& pos) const
{
	// test if pos is inside the cube
	return pos.x >= b[0].x && pos.y >= b[0].y && pos.z >= b[0].z &&
		pos.x <= b[1].x && pos.y <= b[1].y && pos.z <= b[1].z;
}

Scene::Scene()
{
#if TEAPOD
	mesh = new Mesh("assets/teapot.obj");
	bvh = new BVH(mesh);
	return;
#endif
	// setup is done in levels
}


bool Tmpl8::Scene::SphereIntersect(Ray& ray,const Sphere& sphere) const
{
	const float3 sphereDirection = ray.O - sphere.center; // order reversed

	const float a = dot(ray.D, ray.D);
	const float b = 2 * dot(sphereDirection, ray.D);
	const float c = dot(sphereDirection, sphereDirection) - sphere.radius * sphere.radius;

	const float d = b * b - 4 * a * c;
	if (d < 0.0f)
	{
		return false;
	}

	const float discSquared = sqrtf(d);
	const float a2 = 2 * a;
	const float t1 = (-b + discSquared) / a2;
	const float t2 = (-b - discSquared) / a2;

	float dist = std::min(t1, t2);
	if (ray.t > dist && dist > 0) {
		ray.t = dist;
		ray.voxel = sphere.material;
		ray.normal = sphere.center;
		ray.type = Primitive::SpherePrimitive;
		return true;
	}
	return false;
}

bool Tmpl8::Scene::IsSphereOccludeed(const Ray& ray,const Sphere& sphere) const
{
	const float3 sphereDirection = ray.O - sphere.center; // order reversed

	const float a = dot(ray.D, ray.D);
	const float b = 2 * dot(sphereDirection, ray.D);
	const float c = dot(sphereDirection, sphereDirection) - sphere.radius * sphere.radius;

	const float d = b * b - 4 * a * c;
	if (d < 0.0f)
	{
		return false;
	}

	const float discSquared = sqrtf(d);
	const float a2 = 2 * a;
	const float t1 = (-b + discSquared) / a2;
	const float t2 = (-b - discSquared) / a2;

	float dist = std::min(t1, t2);
	if (dist < 0)
	{
		return false;
	}
	return ray.t > dist;
}

int Scene::FindNearest(Ray& ray) const
{
#if TEAPOD
	bvh->Intersect(ray, 0);
	return;
#endif

	int index = -1;
	int spheres = (int)Spheres.size();

	for (int i = 0; i < spheres; i++)
	{
		if (SphereIntersect(ray, Spheres[i]))
		{
			index = i;
		}
	}

	int3 pos;
	int grids = (int)Grids.size();
	for (int i = 0; i < grids; i++)
	{
		if (Grids[i].Traverse(ray, pos))
		{
			index = i;
		}
	}
	
	int planes = (int)Planes.size();
	for (int i = 0; i < planes; i++)
	{
		if (Planes[i].Intersect(ray))
		{
			index = i;
		}
	}

	return index;
}

bool Scene::IsOccluded(Ray& ray) const
{
#if TEAPOD

	// return	bvh->Intersect(ray, 0);

#endif
	
	for (Sphere sphere : Spheres)
	{
		if (IsSphereOccludeed(ray, sphere))
		{
			return true;
		}
	}

	for (Plane plane : Planes)
	{
		plane.Intersect(ray);
	}

	for (VoxelGrid grid : Grids)
	{
		if (grid.IsOcluded(ray))
		{
			return true;
		}
	}

	return false;
}

void Tmpl8::Scene::FindExit(Ray& ray, int index) const
{
	if (index < 0)
	{
		return;
	}

	if (ray.type == Primitive::SpherePrimitive)
	{
		Spheres[index].SphereExit(ray);
	}
	else if (ray.type == Primitive::GridPrimitive)
	{
		Grids[index].FindExit(ray);
	}
}

void Tmpl8::Scene::ClearObjects()
{
	Grids.clear();
	Spheres.clear();
	Planes.clear();
}

VoxelGrid* Tmpl8::Scene::AddVoxelGrid(const VoxelGrid& grid)
{
	Grids.push_back(grid);
	return &Grids.back();
}

Sphere* Tmpl8::Scene::AddSphere(const Sphere& sphere)
{
	Spheres.push_back(sphere);
	return &Spheres.back();
}

Plane* Tmpl8::Scene::AddPlane(const Plane& plane)
{
	Planes.push_back(plane);
	return &Planes.back();
}

VoxelGrid* Tmpl8::Scene::GetVoxelGrid(const int index)
{
	return &Grids[index];
}
