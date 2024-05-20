#pragma once
#include "include/Mesh.h"

struct BVHNode
{
	float3 aabbMin; uint leftFirst;
	float3 aabbMax; uint triCount;
	bool isLeaf() const { return triCount > 0; } // empty BVH leaves do not exist
	float CalculateNodeCost()
	{
		float3 e = aabbMax - aabbMin; // extent of the node
		return (e.x * e.y + e.y * e.z + e.z * e.x) * triCount;
	}
};

class BVH
{
private:
	struct BuildJob
	{
		uint nodeIdx;
		float3 centroidMin, centroidMax;
	};
public:
	BVH() = default;
	BVH(Mesh* mesh);

	void BuildBVH();
	bool Intersect(Ray& ray, uint instanceIdx);
private:
	void Subdivide(uint nodeIdx, uint depth, uint& nodePtr, float3& centroidMin, float3& centroidMax);
	void UpdateNodeBounds(uint nodeIdx, float3& centroidMin, float3& centroidMax);
	float FindBestSplitPlane(BVHNode& node, int& axis, int& splitPos, float3& centroidMin, float3& centroidMax);
	struct Mesh* mesh = 0;
public:
	uint* triIdx = 0;
	uint nodesUsed;
	BVHNode* bvhNode = 0;
	bool subdivToOnePrim = false; // for TLAS experiment
	BuildJob buildStack[64];
	int buildStackPtr;
};