#pragma once

struct SpotLight
{
	float4 dummy; // 16 bytes
	float3 position; // 12 bytes
	float3 direction; // 12 bytes
	float3 color; // 12 bytes
	float angle; // 4 bytes
};