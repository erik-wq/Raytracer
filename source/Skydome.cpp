#include "precomp.h"
#include "include/Skydome.h"
#include <stb_image.h>

Skydome::Skydome()
{
	int chanels;
	pixels = stbi_loadf("assets/lilienstein_4k.hdr", &width, &height, &chanels, 0);
	for (int i = 0; i < width * height * 3; i++)
	{
		pixels[i] = sqrtf(pixels[i]);
	}
}

Skydome::~Skydome()
{
	stbi_image_free(pixels);
}

float3 Skydome::GetSkyColor(Ray& ray)
{
	float3 direction = ray.D;

	// Sample sky
	const float uFloat = static_cast<float>(width) * atan2f(direction.z, direction.x) * INV2PI - 0.5f;
	const float vFloat = static_cast<float>(height) * acosf(direction.y) * INVPI - 0.5f;

	const int u = static_cast<int>(uFloat);
	const int v = static_cast<int>(vFloat);

	const int skyIdx = max(0, u + v * width);
	return float3(pixels[skyIdx * 3], pixels[skyIdx * 3 + 1], pixels[skyIdx * 3 + 2]);
}
