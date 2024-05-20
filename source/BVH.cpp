#include "precomp.h"
#include "include/BHV.h"

BVH::BVH(Mesh* mesh)
{
	this->mesh = mesh;
	bvhNode = new BVHNode[2 * mesh->numTriangles - 1];
	triIdx = new uint[mesh->numTriangles];
	BuildBVH();
}

void BVH::BuildBVH()
{
	nodesUsed = 2;
	const int count = mesh->numTriangles;

	memset(bvhNode, 0, (count * 2 - 1) * sizeof(BVHNode));

	// init triangle indices
	for (int i = 0; i < count; i++)
	{
		triIdx[i] = i;
	}

	Triangle* triangle = mesh->triangles;
	for (int i = 0; i < count; i++)
	{
		// calculating centroid of the triangle
		triangle[i].controid = (triangle[i].v0 + triangle[i].v1 + triangle[i].v2) * 0.3333f;
	}


	BVHNode& root = bvhNode[0];
	// init root node
	root.leftFirst = 0, root.triCount = count;
	float3 centroidMin, centroidMax;
	UpdateNodeBounds(0, centroidMin, centroidMax);

	// recursive subdivision
	buildStackPtr = 0;
	Subdivide(0, 0, nodesUsed, centroidMin, centroidMax);
}

bool BVH::Intersect(Ray& ray, uint)
{
	BVHNode* node = &bvhNode[0], * stack[64];
	uint stackPtr = 0;
	bool intersected = false;
	while (1)
	{
		if (node->isLeaf())
		{
			float3 normal = {0};
			for (uint i = 0; i < node->triCount; i++)
			{
				uint instPrim = triIdx[node->leftFirst + i];
				// final triangle mesh intersection
				intersected |= IntersectTriangle(ray, mesh->triangles[instPrim], mesh->normals[instPrim]);
			}

			if (stackPtr == 0)
			{
				break;
			}
			else
			{
				node = stack[--stackPtr];
				continue;
			}
		}
		BVHNode* child1 = &bvhNode[node->leftFirst];
		BVHNode* child2 = &bvhNode[node->leftFirst + 1];

		float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax);
		float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax);

		if (dist1 > dist2) 
		{ 
			swap(dist1, dist2); 
			swap(child1, child2); 
		}
		if (dist1 == 1e30f)
		{
			if (stackPtr == 0) break; else node = stack[--stackPtr];
		}
		else
		{
			node = child1;
			if (dist2 != 1e30f) stack[stackPtr++] = child2;
		}
	}
	return intersected;
}

void BVH::Subdivide(uint nodeIdx, uint depth, uint& nodePtr, float3& centroidMin, float3& centroidMax)
{
	BVHNode& node = bvhNode[nodeIdx];
	// determine split axis using SAH
	int axis, splitPos;
	float splitCost = FindBestSplitPlane(node, axis, splitPos, centroidMin, centroidMax);
	// terminate recursion
	if (subdivToOnePrim)
	{
		if (node.triCount == 1) return;
	}
	else
	{
		float nosplitCost = node.CalculateNodeCost();
		if (splitCost >= nosplitCost) return;
	}
	// in-place partition
	int i = node.leftFirst;
	int j = i + node.triCount - 1;
	float scale = 8 / (centroidMax[axis] - centroidMin[axis]);
	while (i <= j)
	{
		// use the exact calculation we used for binning to prevent rare inaccuracies
		int binIdx = min(8 - 1, (int)((mesh->triangles[triIdx[i]].controid[axis] - centroidMin[axis]) * scale));
		if (binIdx < splitPos) i++; else swap(triIdx[i], triIdx[j--]);
	}
	// abort split if one of the sides is empty
	int leftCount = i - node.leftFirst;
	if (leftCount == 0 || leftCount == (int)node.triCount) return; // never happens for dragon mesh, nice
	// create child nodes
	int leftChildIdx = nodePtr++;
	int rightChildIdx = nodePtr++;
	bvhNode[leftChildIdx].leftFirst = node.leftFirst;
	bvhNode[leftChildIdx].triCount = leftCount;
	bvhNode[rightChildIdx].leftFirst = i;
	bvhNode[rightChildIdx].triCount = node.triCount - leftCount;
	node.leftFirst = leftChildIdx;
	node.triCount = 0;
	// recurse
	UpdateNodeBounds(leftChildIdx, centroidMin, centroidMax);
	Subdivide(leftChildIdx, depth + 1, nodePtr, centroidMin, centroidMax);
	UpdateNodeBounds(rightChildIdx, centroidMin, centroidMax);
	Subdivide(rightChildIdx, depth + 1, nodePtr, centroidMin, centroidMax);
}

void BVH::UpdateNodeBounds(uint nodeIdx, float3& centroidMin, float3& centroidMax)
{
	// node reference
	BVHNode& node = bvhNode[nodeIdx];

	node.aabbMin = float3(1e30f);
	node.aabbMax = float3(-1e30f);
	centroidMin = float3(1e30f);
	centroidMax = float3(-1e30f);
	for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
	{
		uint leafTriIdx = triIdx[first + i];
		Triangle& leafTri = mesh->triangles[leafTriIdx];
		node.aabbMin = fminf(node.aabbMin, leafTri.v0);
		node.aabbMin = fminf(node.aabbMin, leafTri.v1);
		node.aabbMin = fminf(node.aabbMin, leafTri.v2);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.v0);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.v1);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.v2);
		centroidMin = fminf(centroidMin, leafTri.controid);
		centroidMax = fmaxf(centroidMax, leafTri.controid);
	}
}

float BVH::FindBestSplitPlane(BVHNode& node, int& axis, int& splitPos, float3& centroidMin, float3& centroidMax)
{
	// starting cost
	float bestCost = 1e30f;

	for (int a = 0; a < 3; a++)
	{
		float boundsMin = centroidMin[a], boundsMax = centroidMax[a];
		if (boundsMin == boundsMax) continue;

		// populate the bins BINS = 8
		float scale = 8 / (boundsMax - boundsMin);
		float leftCountArea[8 - 1], rightCountArea[8 - 1];
		int leftSum = 0, rightSum = 0;

		struct Bin { aabb bounds; int triCount = 0; } bin[8];
		for (uint i = 0; i < node.triCount; i++)
		{
			Triangle& triangle = mesh->triangles[triIdx[node.leftFirst + i]];
			int binIdx = min(8 - 1, (int)((triangle.controid[a] - boundsMin) * scale));
			bin[binIdx].triCount++;
			bin[binIdx].bounds.Grow(triangle.v0);
			bin[binIdx].bounds.Grow(triangle.v1);
			bin[binIdx].bounds.Grow(triangle.v2);
		}
		// gather data for the 7 planes between the 8 bins
		aabb leftBox, rightBox;
		for (int i = 0; i < 8 - 1; i++)
		{
			leftSum += bin[i].triCount;
			leftBox.Grow(bin[i].bounds);
			leftCountArea[i] = leftSum * leftBox.Area();
			rightSum += bin[8 - 1 - i].triCount;
			rightBox.Grow(bin[8 - 1 - i].bounds);
			rightCountArea[8 - 2 - i] = rightSum * rightBox.Area();
		}

		// calculate SAH cost for the 7 planes
		scale = (boundsMax - boundsMin) / 8;
		for (int i = 0; i < 8 - 1; i++)
		{
			const float planeCost = leftCountArea[i] + rightCountArea[i];
			if (planeCost < bestCost)
				axis = a, splitPos = i + 1, bestCost = planeCost;
		}
	}

	return bestCost;
}
