#pragma once

struct DirectionalLight
{
public:
	float3 direction = float3(0, -1, 0); // 12 bytes
	float3 color = float3(1); // 12 bytes
	double dummy; // 8 bytes
};