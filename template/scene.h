#pragma once

// high level settings
// #define TWOLEVEL
#define WORLDSIZE 256// power of 2. Warning: max 512 for a 512x512x512x4 bytes = 512MB world!
// #define USE_SIMD
// #define USE_FMA3
// #define SKYDOME
// #define WHITTED
// #define DOF

// low-level / derived
#define WORLDSIZE2	(WORLDSIZE*WORLDSIZE)
#ifdef TWOLEVEL
// feel free to replace with whatever suits your two-level implementation,
// should you chose this challenge.
#define BRICKSIZE	8
#define BRICKSIZE2	(BRICKSIZE*BRICKSIZE)
#define BRICKSIZE3	(BRICKSIZE*BRICKSIZE*BRICKSIZE)
#define GRIDSIZE	(WORLDSIZE/BRICKSIZE)
#define VOXELSIZE	(1.0f/WORLDSIZE)
#else
#define GRIDSIZE	WORLDSIZE
#endif
#define GRIDSIZE2	(GRIDSIZE*GRIDSIZE)
#define GRIDSIZE3	(GRIDSIZE*GRIDSIZE*GRIDSIZE)

#define TEAPOD 0

class BVH;
class VoxelGrid;
class Sphere;
class Plane;

namespace Tmpl8 {

class Ray
{
public:
	Ray() = default;
	Ray( const float3 origin, const float3 direction, const float rayLength = 1e34f, const int rgb = 0 )
		: O( origin ), D( direction ), t( rayLength ), voxel( rgb )
	{
		if (D.x != 0)
		{
			rD.x = 1 / D.x;
		}
		if (D.x != 0)
		{
			rD.y = 1 / D.y;
		}
		if (D.z != 0)
		{
			rD.z = 1 / D.z;
		}
	}
	float3 IntersectionPoint() const { return O + t * D; }
	float3 GetNormal() const;
	// ray data
	float3 O;					// ray origin
	float3 rD = {0};					// reciprocal ray direction
	float3 D = {0};		// ray direction
	float t = 1e34f;			// ray length
	// float3 Dsign = float3( 1 );	// inverted ray direction signs, -1 or 1
	uint voxel = 0;			// 8 bit value, used for voxel coloring
	float3 normal = float3( 0 );	// normal at intersection point
	Primitive type;				// type of intersected primitive

	__inline static float3 min3(const float3& a, const float3& b)
	{
		return float3(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
	}
};

// intersect methods
inline float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax)
{
	float tx1 = (bmin.x - ray.O.x) * ray.rD.x, tx2 = (bmax.x - ray.O.x) * ray.rD.x;
	float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
	float ty1 = (bmin.y - ray.O.y) * ray.rD.y, ty2 = (bmax.y - ray.O.y) * ray.rD.y;
	tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
	float tz1 = (bmin.z - ray.O.z) * ray.rD.z, tz2 = (bmax.z - ray.O.z) * ray.rD.z;
	tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
	if (tmax >= tmin && tmin < ray.t && tmax > 0) return tmin; else return 1e30f;
};

inline bool IntersectTriangle(Ray& ray, const Triangle& triangle, const TriangleNormal&)
{
	const float3 e1 = triangle.v1 - triangle.v0;
	const float3 e2 = triangle.v2 - triangle.v0;
	const float3 h = cross(ray.D, e2);
	const float a = dot(e1, h);
	if (a > -0.00001f  && a < 0.00001f)
	{
		return false;
	}
	const float f = 1.0f / a;
	const float3 s = ray.O - triangle.v0;
	float u = f * dot(s, h);
	if (u < 0.0f || u > 1.0f)
	{
		return false;
	}
	const float3 q = cross(s, e1);
	float v = f * dot(ray.D, q);
	if (v < 0.0f || u + v > 1.0f)
	{
		return false;
	}
	const float t = f * dot(e2, q);
	if (ray.t > t)
	{
		ray.t = t;
		ray.voxel = 1;
		// ray.normal = u * normal.n1 + v * normal.n2 + (1 - (u + v)) * normal.n0;
	}
	return true;
}

class Cube
{
public:
	Cube() = default;
	Cube( const float3 pos, const float3 size );
	float Intersect( const Ray& ray ) const;
	bool Contains( const float3& pos ) const;
	float3 b[2];
};

class Scene
{
public:
	struct DDAState
	{
		int3 step;				// 16 bytes
		uint X, Y, Z;			// 12 bytes
		float t;				// 4 bytes
		float3 tdelta;
		float dummy1 = 0;		// 16 bytes
		float3 tmax;
		float dummy2 = 0;		// 16 bytes, 64 bytes in total
	};
	Scene();
	int FindNearest( Ray& ray ) const;
	bool IsOccluded( Ray& ray ) const;
	void FindExit(Ray& ray, int index) const;

	void ClearObjects();

	VoxelGrid* AddVoxelGrid(const VoxelGrid& grid);
	Sphere* AddSphere(const Sphere& sphere);
	Plane* AddPlane(const Plane& plane);

	VoxelGrid* GetVoxelGrid(const int index);

private:

	std::vector<VoxelGrid> Grids;
	std::vector<Sphere> Spheres;
	std::vector<Plane> Planes;

	bool SphereIntersect( Ray& ray,const Sphere& sphere) const;
	bool IsSphereOccludeed(const Ray& ray,const Sphere& sphere ) const;
};

}