#include "precomp.h"
#include "include/Mesh.h"
#include "lib/tiny_obj_loader.h"

Mesh::Mesh(const char* filename)
{
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename, "");

	if (!ret)
	{
		return;
	}

	std::vector<Triangle> trian;
	std::vector<Triangle> finalTriangles;


	std::vector<TriangleNormal> norm;
	std::vector<TriangleNormal> finalNorm;

	// can add nromals and texture coordinates
	for (const auto shape : shapes)
	{
		trian.resize(trian.size() + shape.mesh.indices.size() / 3);
		norm.resize(trian.size() + shape.mesh.indices.size() / 3);
		int count = 0;
		for (size_t f = 0; f < shape.mesh.indices.size() / 3; f++) {
			tinyobj::index_t idx0 = shape.mesh.indices[3 * f + 0];
			tinyobj::index_t idx1 = shape.mesh.indices[3 * f + 1];
			tinyobj::index_t idx2 = shape.mesh.indices[3 * f + 2];

			Triangle* curr = &trian[count];
			int f0 = idx0.vertex_index;
			int f1 = idx1.vertex_index;
			int f2 = idx2.vertex_index;

			for (int k = 0; k < 3; k++) {
				curr->v0.cell[k] = attrib.vertices[3 * f0 + k];
				curr->v1.cell[k] = attrib.vertices[3 * f1 + k];
				curr->v2.cell[k] = attrib.vertices[3 * f2 + k];
			}

			f0 = idx0.normal_index;
			f1 = idx1.normal_index;
			f2 = idx2.normal_index;

			TriangleNormal* currNorm = &norm[count];
			if (attrib.normals.size() > 0) {
				for (int k = 0; k < 3; k++) {
					currNorm->n0.cell[k] = attrib.normals[3 * f0 + k];
					currNorm->n1.cell[k] = attrib.normals[3 * f1 + k];
					currNorm->n2.cell[k] = attrib.normals[3 * f2 + k];
				}
			}

			count++;
			numTriangles++;
		}
		finalTriangles.insert(finalTriangles.end(), trian.begin(), trian.end());
		finalNorm.insert(finalNorm.end(), norm.begin(), norm.end());
	}

	triangles = new Triangle[finalTriangles.size()];
	normals = new TriangleNormal[finalNorm.size()];

	memcpy(triangles, finalTriangles.data(), finalTriangles.size() * sizeof(Triangle));
	memcpy(normals, finalNorm.data(), finalNorm.size() * sizeof(TriangleNormal));
}
