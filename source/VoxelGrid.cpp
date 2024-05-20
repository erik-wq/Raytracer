#include "precomp.h"
#include "include/VoxelGrid.h"
#include <unordered_map>

#define OGT_VOX_IMPLEMENTATION
#include "lib/ogt_vox.h"

VoxelGrid::VoxelGrid()
{
	// the voxel world sits in a 1x1x1 cube
	cube = Cube(float3(0, 0, 0), float3(1, 1, 1));
	// initialize the scene using Perlin noise, parallel over z

	for (int i = 0; i < 4; i++)
	{
		cellSizes[i] = 1.0f / (float)sizes[i];
	}

	grid = (uint8_t*)MALLOC64(TOTALCELLS * sizeof(uint8_t));
	memset(grid, 0, TOTALCELLS * sizeof(uint8_t));

#if 1
	ReadVoxFile();
	return;
#else

	if (ReadGridData())
	{
		return;
	}
	
	GenerateGridData();
#endif
}
void VoxelGrid::RemoveVoxel(const int3& pos)
{
	int offset = offsets[0];
	int size = sizes[0];
	int size2 = sizes2[0];

	// set lod0 to 0
	grid[offset + pos.x + pos.y * size + pos.z * size2] = 0;
	
	// check higher lod
	int3 cellPos;

	// Check Higher LODs
	for (int i = 1; i < 4; i++)
	{
		offset = offsets[i-1];
		size = sizes[i-1];
		size2 = sizes2[i-1];

		cellPos.x = (pos.x >> i) << 1;
		cellPos.y = (pos.y >> i) << 1;
		cellPos.z = (pos.z >> i) << 1;

		for (int x = 0; x < 2; x++)
		{
			for (int y = 0; y < 2; y++)
			{
				for (int z = 0; z < 2; z++)
				{
					if (grid[offset + cellPos.x + x + (cellPos.y + y) * size + (cellPos.z + z) * size2])
					{
						return;
					}
				}
			}
		}

		offset = offsets[i];
		size = sizes[i];
		size2 = sizes2[i];
		grid[offset + (pos.x >> i) + (pos.y >> i)* size + (pos.z >> i) * size2] = 0;
	}
}
void VoxelGrid::ReadVoxFile()
{
	// load via opengametools .vox file reader
	FILE* fp = fopen("./assets/maze.vox", "rb");
	uint32_t buffer_size = _filelength(_fileno(fp));
	uint8_t* buffer = new uint8_t[buffer_size];
	fread(buffer, buffer_size, 1, fp);
	fclose(fp);
	const ogt_vox_scene* scene = ogt_vox_read_scene(buffer, buffer_size);
	delete[] buffer;
	// find bounds of voxel scene
	int3 bmin = make_int3(0xffffff), bmax = make_int3(-0xffffff);
	for (uint i = 0; i < scene->num_instances; i++)
	{
		const ogt_vox_instance& instance = scene->instances[i];
		const mat4 M = mat4::FromColumnMajor(*(mat4*)&instance.transform);
		int3 pos = make_int3(M.GetTranslation());
		const ogt_vox_model* model = scene->models[instance.model_index];
		int3 size = make_int3(M * make_float4(make_int3(model->size_x, model->size_y, model->size_z), 0));
		bmin = min(bmin, min(pos, pos + size));
		bmax = max(bmax, max(pos, pos + size));
	}

	std::unordered_map<uint8_t, uint8_t> data;
	// store the data in the frame
	for (uint n = 0; n < scene->num_instances; n++)
	{
		// get instance model data
		// note: magicavoxel author is evil; voxel data may be transformed.
		const ogt_vox_instance& instance = scene->instances[n];
		const mat4 M = mat4::FromColumnMajor(*(mat4*)&instance.transform);
		const int3 pos = make_int3(M.GetTranslation()) - bmin;
		const ogt_vox_model* model = scene->models[instance.model_index];
		const int3 size = make_int3(model->size_x, model->size_y, model->size_z);
		const int3 MX = make_int3(make_float3(M.cell[0], M.cell[1], M.cell[2]));
		const int3 MY = make_int3(make_float3(M.cell[4], M.cell[5], M.cell[6]));
		const int3 MZ = make_int3(make_float3(M.cell[8], M.cell[9], M.cell[10]));

		// store actual voxels
		for (int i = 0, z = 0; z < size.z; z++)
			for (int y = 0; y < size.y; y++)
				for (int x = 0; x < size.x; x++, i++)
				{
					const uint8_t v = model->voxel_data[i];
					if (!v) continue;
					uint8_t mat = 255;

					if (v == 211) mat = 1;
					else if (v == 226) mat = 2;
					else if (v == 112) mat = 3;
					else if (v == 35) mat = 5;

					if (mat == 255) continue;

					SetLOD(x, z, y, mat);
				}
	}

	ogt_vox_destroy_scene(scene);
}


void VoxelGrid::SetLOD(const uint x, const uint y, const uint z, const uint8_t v)
{
	int pos = 0;

	//LOD1
	pos = x / 2 + y / 2 * LOD1SIZE + z / 2 * LOD12;
	pos += LOD1OFFSET;
	grid[pos] |= v;

	//LOD2
	pos = x / 4 + y / 4 * LOD2SIZE + z / 4 * LOD22;
	pos += LOD2OFFSET;
	grid[pos] |= v;

	// LOD3
	pos = x / 8 + y / 8 * LOD3SIZE + z / 8 * LOD32;
	grid[pos] |= v;

	// LOD0
	pos = GRIDOFFSET + x + y * GRIDSIZE + z * GRIDSIZE2;
	grid[pos] = v;
}

bool VoxelGrid::ReadGridData()
{
	FILE* file = fopen("scene.save", "rb");
	if (!file)
	{
		return false;
	}

	fseek(file, 0, SEEK_END);
	long file_size = ftell(file);
	rewind(file);

	size_t result = fread(grid, 1, file_size, file);
	if (result != file_size) {
		fclose(file);
		return false;
	}

	fclose(file);
	return true;
}

void VoxelGrid::GenerateGridData()
{
	for (int z = 0; z < WORLDSIZE; z++)
	{
		const float fz = (float)z / WORLDSIZE;
		for (int y = 0; y < WORLDSIZE; y++)
		{
			const float fy = (float)y / WORLDSIZE;
			float fx = 0;
			for (int x = 0; x < WORLDSIZE; x++, fx += 1.0f / WORLDSIZE)
			{
				const float n = noise3D(fx, fy, fz);
				uint8_t v = (uint8_t)Rand(MaxMaterail) + 1;
				SetLOD(x, y, z, n > 0.09f ? v : 0x00);
			}
		}
	}

	FILE* file = fopen("scene.save", "wb");
	fwrite(grid, 1, TOTALCELLS, file);
	fclose(file);
}

bool VoxelGrid::Setup3DDDA(const Ray& ray, DDAState& state, int lod, const float3& sign) const
{

	// if ray is not inside the world: advance until it is
	state.t = 0;
	if (!cube.Contains(ray.O))
	{
		state.t = cube.Intersect(ray);
		if (state.t > 1e33f) return false; // ray misses voxel data entirely
	}

	float size = (float)sizes[lod];
	state.gridSize = (uint)size;

	float cellSize = cellSizes[lod];
	state.step = make_int3(1 - sign * 2);
	const float3 posInGrid = (float)size * (ray.O + (state.t + 0.00005f) * ray.D);
	const float3 gridPlanes = (ceilf(posInGrid) - sign) * cellSize;
	const int3 P = clamp(make_int3(posInGrid), 0, sizes[lod] - 1);
	state.X = P.x, state.Y = P.y, state.Z = P.z;
	state.tdelta = cellSize * float3(state.step) * ray.rD;
	state.tmax = (gridPlanes - ray.O) * ray.rD;
	// proceed with traversal

	return true;
}


bool VoxelGrid::SSESetup3DDDA(const SSERay& ray, const Ray& normalRay , DDAState& state, int lod) const
{
	// if ray is not inside the world: advance until it is
	state.t = 0;
	if (!cube.Contains(normalRay.O))
	{
		state.t = cube.Intersect(normalRay);
		if (state.t > 1e33f) return false; // ray mis1ses voxel data entirely
	}

	float size = (float)sizes[lod];
	state.gridSize = (uint)size;

	__m128 t = _mm_set1_ps(state.t + 0.00005f);
	__m128 gridSize = _mm_set1_ps(size);

	__m128 cell = _mm_set1_ps(cellSizes[lod]);

	__m128 step4 = _mm_sub_ps(one, _mm_mul_ps(ray.sign, two));
	__m128 posInGrid4 = _mm_mul_ps(gridSize, _mm_add_ps(ray.o, _mm_mul_ps(t, ray.d)));
	__m128 P4 = _mm_max_ps(_mm_min_ps(posInGrid4, _mm_sub_ps(gridSize, one)), zero);
	__m128 gridPlane = _mm_mul_ps(_mm_set1_ps(cellSizes[lod]), _mm_sub_ps(_mm_ceil_ps(posInGrid4), ray.sign));

	// __m128 rD4 = _mm_div_ps(one, D4);  //_mm_rcp_ps(D4); aprochimate reciprical isn't enought
	float4 delta4, tmax4;
	_mm_store_ps(delta4.cell, _mm_mul_ps(cell, _mm_mul_ps(step4, ray.rD)));
	_mm_store_ps(tmax4.cell, _mm_mul_ps(_mm_sub_ps(gridPlane, ray.o), ray.rD));

	state.step = make_int3((int)step4.m128_f32[0], (int)step4.m128_f32[1], (int)step4.m128_f32[2]);
	const int3 P = make_int3((int)P4.m128_f32[0], (int)P4.m128_f32[1], (int)P4.m128_f32[2]);
	state.X = P.x, state.Y = P.y, state.Z = P.z;
	state.tdelta = (float3)delta4;
	state.tmax = (float3)tmax4;
	return true;
}


void VoxelGrid::ChangeLOD(Ray& ray, DDAState& state, int lod, const float3& sign) const
{
	float size = (float)sizes[lod];
	state.gridSize = (uint)size;

	float cellSize = cellSizes[lod];

	const float3 posInGrid = (float)size * (ray.O + (state.t + 0.00005f) * ray.D);
	const float3 gridPlanes = (ceilf(posInGrid) - sign) * cellSize;
	const int3 P = clamp(make_int3(posInGrid), 0, sizes[lod] - 1);
	state.X = P.x, state.Y = P.y, state.Z = P.z;
	state.tdelta *= 0.5f;
	state.tmax = (gridPlanes - ray.O) * ray.rD;
}

void VoxelGrid::SSEChangeLOD(const SSERay& ray, DDAState& state, int lod) const
{
	float size = (float)sizes[lod];
	state.gridSize = (uint)size;

	__m128 t = _mm_set1_ps(state.t + 0.00005f);
	__m128 gridSize = _mm_set1_ps(size);

	__m128 cell = _mm_set1_ps(cellSizes[lod]);

	__m128 posInGrid4 = _mm_mul_ps(gridSize, _mm_add_ps(ray.o, _mm_mul_ps(t, ray.d)));
	__m128 P4 = _mm_max_ps(_mm_min_ps(posInGrid4, _mm_sub_ps(gridSize, one)), zero);
	__m128 gridPlane = _mm_mul_ps(_mm_set1_ps(cellSizes[lod]), _mm_sub_ps(_mm_ceil_ps(posInGrid4), ray.sign));

	// __m128 rD4 = _mm_div_ps(one, D4);  //_mm_rcp_ps(D4); aprochimate reciprical isn't enought
	float4 tmax4;
	_mm_store_ps(tmax4.cell, _mm_mul_ps(_mm_sub_ps(gridPlane, ray.o), ray.rD));

	const int3 P = make_int3((int)P4.m128_f32[0], (int)P4.m128_f32[1], (int)P4.m128_f32[2]);
	state.X = P.x, state.Y = P.y, state.Z = P.z;
	state.tdelta *= 0.5f;
	state.tmax = (float3)tmax4;
}

void VoxelGrid::FindExit(Ray& ray) const
{
	// cannont use multi level grid since higher cell could skip some empty cells
	uint x = (*(uint*)&ray.D.x >> 31);
	uint y = (*(uint*)&ray.D.y >> 31);
	uint z = (*(uint*)&ray.D.z >> 31);

	float3 Dsign = float3((float)x, (float)y, (float)z);

	DDAState s;
	if (!Setup3DDDA(ray, s, 3, Dsign)) return;
	// start stepping
	while (1)
	{
		const uint cell = grid[s.X + s.Y * GRIDSIZE + s.Z * GRIDSIZE2];
		if (!cell)
		{
			ray.t = s.t;
			ray.type = Primitive::GridPrimitive;
			break;
		}
		if (s.tmax.x < s.tmax.y)
		{
			if (s.tmax.x < s.tmax.z)
			{
				s.t = s.tmax.x, s.X += s.step.x;
				if (s.X >= GRIDSIZE)
					break;
				s.tmax.x += s.tdelta.x;
			}
			else
			{
				s.t = s.tmax.z, s.Z += s.step.z;
				if (s.Z >= GRIDSIZE)
				{
					break;
				}
				s.tmax.z += s.tdelta.z;
			}
		}
		else
		{
			if (s.tmax.y < s.tmax.z)
			{
				s.t = s.tmax.y, s.Y += s.step.y;
				if (s.Y >= GRIDSIZE)
				{

					break;
				}
				s.tmax.y += s.tdelta.y;
			}
			else
			{
				s.t = s.tmax.z, s.Z += s.step.z;
				if (s.Z >= GRIDSIZE)
				{
					break;
				}
				s.tmax.z += s.tdelta.z;
			}
		}
	}
}


// 512 total L1 cells
bool VoxelGrid::Traverse(Ray& ray, int3& outPos) const
{
	uint x = (*(uint*)&ray.D.x >> 31);
	uint y = (*(uint*)&ray.D.y >> 31);
	uint z = (*(uint*)&ray.D.z >> 31);

	float3 Dsign = float3((float)x, (float)y, (float)z);

#if SSE
	SSERay sseRay;
	sseRay.o = _mm_set_ps(1, ray.O.z, ray.O.y, ray.O.x);
	sseRay.d = _mm_set_ps(1, ray.D.z, ray.D.y, ray.D.x);
	sseRay.rD = _mm_div_ps(one, sseRay.d);
	sseRay.sign = _mm_set_ps(1, Dsign.z, Dsign.y, Dsign.x);
#endif

	int lod = 3;
	DDAState s;
#if  SSE
	if (!SSESetup3DDDA(sseRay, ray, s, lod)) return false;
#else
	if (!Setup3DDDA(ray, s, lod, Dsign)) return false;
#endif
	int offset = offsets[lod];
	int size = sizes[lod];
	int size2 = sizes2[lod];

	// start stepping
	while (s.t < ray.t)
	{
		const uint cell = grid[offset + s.X + s.Y * size + s.Z * size2];
		if (cell)
		{
			if (lod == 0)
			{
				if (ray.t < s.t)
				{
					return false;
				}

				ray.t = s.t;
				ray.voxel = cell;
				ray.type = Primitive::GridPrimitive;
				outPos = make_int3(s.X, s.Y, s.Z);
				return true;
			}
			else
			{
				lod--;
				offset = offsets[lod];
				size = sizes[lod];
				size2 = sizes2[lod];
#if  SSE
				SSEChangeLOD(sseRay, s, lod);
#else
				ChangeLOD(ray, s, lod, Dsign);
#endif
				continue;
			}
		}

		if (s.tmax.x < s.tmax.y)
		{
			if (s.tmax.x < s.tmax.z)
			{
				if ((s.X += s.step.x) >= s.gridSize)
					return false;
				s.t = s.tmax.x;
				s.tmax.x += s.tdelta.x;
			}
			else
			{
				if ((s.Z += s.step.z) >= s.gridSize)
					return false;

				s.t = s.tmax.z;
				s.tmax.z += s.tdelta.z;
			}
		}
		else
		{
			if (s.tmax.y < s.tmax.z)
			{
				if ((s.Y += s.step.y) >= s.gridSize)
					return false;

				s.t = s.tmax.y;
				s.tmax.y += s.tdelta.y;
			}
			else
			{
				if ((s.Z += s.step.z) >= s.gridSize)
					return false;

				s.t = s.tmax.z;
				s.tmax.z += s.tdelta.z;
			}
		}
	}
	return false;
}


bool VoxelGrid::IsOcluded(Ray& ray) const
{
	uint x = (*(uint*)&ray.D.x >> 31);
	uint y = (*(uint*)&ray.D.y >> 31);
	uint z = (*(uint*)&ray.D.z >> 31);

	float3 Dsign = float3((float)x, (float)y, (float)z);

#if SSE
	SSERay sseRay;
	sseRay.o = _mm_set_ps(1, ray.O.z, ray.O.y, ray.O.x);
	sseRay.d = _mm_set_ps(1, ray.D.z, ray.D.y, ray.D.x);
	sseRay.rD = _mm_div_ps(one, sseRay.d);
	sseRay.sign = _mm_set_ps(1, Dsign.z, Dsign.y, Dsign.x);
#endif

	int lod = 3;
	DDAState s;
#if  SSE
	if (!SSESetup3DDDA(sseRay, ray, s, lod)) return false;
#else
	if (!Setup3DDDA(ray, s, lod, Dsign)) return false;
#endif
	int offset = offsets[lod];
	int size = sizes[lod];
	int size2 = sizes2[lod];

	// start stepping
	while (s.t < ray.t)
	{
		const uint cell = grid[offset + s.X + s.Y * size + s.Z * size2];
		if (cell)
		{
			if (lod == 0)
			{
				return ray.t > s.t;
			}
			else
			{
				lod--;
				offset = offsets[lod];
				size = sizes[lod];
				size2 = sizes2[lod];
#if  SSE
				SSEChangeLOD(sseRay, s, lod);
#else
				ChangeLOD(ray, s, lod, Dsign);
#endif
				continue;
			}
		}

		if (s.tmax.x < s.tmax.y)
		{
			if (s.tmax.x < s.tmax.z)
			{
				if ((s.X += s.step.x) >= s.gridSize)
					return false;
				s.t = s.tmax.x;
				s.tmax.x += s.tdelta.x;
			}
			else
			{
				if ((s.Z += s.step.z) >= s.gridSize)
					return false;

				s.t = s.tmax.z;
				s.tmax.z += s.tdelta.z;
			}
		}
		else
		{
			if (s.tmax.y < s.tmax.z)
			{
				if ((s.Y += s.step.y) >= s.gridSize)
					return false;

				s.t = s.tmax.y;
				s.tmax.y += s.tdelta.y;
			}
			else
			{
				if ((s.Z += s.step.z) >= s.gridSize)
					return false;

				s.t = s.tmax.z;
				s.tmax.z += s.tdelta.z;
			}
		}
	}
	return false;
}