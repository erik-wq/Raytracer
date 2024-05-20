#pragma once

struct PointLight
{
	float3 position = float3(0); // 12 bytes
	float3 color = float3(0); // 12 bytes
	double dummy = 0; // 8 bytes
};