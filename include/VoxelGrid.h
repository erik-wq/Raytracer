#pragma once

// helper macros for lod
// lod level merges 4 voxels into one
#define LOD3SIZE GRIDSIZE / 8
#define LOD2SIZE GRIDSIZE / 4
#define LOD1SIZE GRIDSIZE / 2

#define LOD32 LOD3SIZE * LOD3SIZE
#define LOD22 LOD2SIZE * LOD2SIZE
#define LOD12 LOD1SIZE * LOD1SIZE

#define LOD3COUNT LOD32 * LOD3SIZE
#define LOD2COUNT LOD22 * LOD2SIZE
#define LOD1COUNT LOD12 * LOD1SIZE

#define LOD3CELLSIZE 1.f / LOD3SIZE
#define LOD2CELLSIZE 1.f / LOD2SIZE
#define LOD1CELLSIZE 1.f / LOD1SIZE

#define GRIDOFFSET LOD3COUNT + LOD2COUNT + LOD1COUNT
#define LOD1OFFSET LOD3COUNT + LOD2COUNT
#define LOD2OFFSET LOD3COUNT

#define TOTALCELLS GRIDSIZE3 + LOD1COUNT + LOD2COUNT + LOD3COUNT

#define SSE 1

struct DDAState
{
	int3 step;				// 16 bytes
	uint X, Y, Z;			// 12 bytes
	float t;				// 4 bytes
	float3 tdelta;
	int lod = 3;
	float3 tmax;
	uint gridSize;
};

struct SSERay
{
	__m128 o;
	__m128 d;
	__m128 rD;
	__m128 sign;
};

class VoxelGrid
{
public:
	VoxelGrid();

	bool IsOcluded(Ray& ray) const;

	bool Traverse(Ray& ray, int3& outPos) const;

	void FindExit(Ray& ray) const;

	void RemoveVoxel(const int3& pos);
private:

	void ReadVoxFile();
	void SetLOD(const uint x, const uint y, const uint z, const uint8_t v);

	bool ReadGridData();

	void GenerateGridData();

	bool Setup3DDDA(const Ray& ray, DDAState& state, int lod, const float3& sign) const;

	void ChangeLOD(Ray& ray, DDAState& state, int lod, const float3& sign) const;

	bool SSESetup3DDDA(const SSERay& ray, const Ray& normalRay, DDAState& state, int lod) const;

	void SSEChangeLOD(const SSERay& ray, DDAState& state, int lod) const;

	__inline static float3 min3(const float3& a, const float3& b)
	{
		return float3(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
	}

	Cube cube;

	uint8_t* grid; 

	__m128 zero = _mm_set1_ps(0);
	__m128 one = _mm_set1_ps(1.0f);
	__m128 two = _mm_set1_ps(2.f);

	int offsets[4] = { GRIDOFFSET, LOD1OFFSET, LOD2OFFSET, 0};
	int sizes[4] = { GRIDSIZE, LOD1SIZE, LOD2SIZE, LOD3SIZE };
	int sizes2[4] = { GRIDSIZE2, LOD12, LOD22, LOD32};
	float cellSizes[4] = {};
};