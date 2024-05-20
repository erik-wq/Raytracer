#pragma once

struct Materail
{
public:
	float4 dummy; 
	float3 albedo = float3(0.5f); 
	float3 AbsorbedColor = float3(0.65f, 0, 0);
	MaterialType type = MaterialType::Diffuse;
	float reflectance = 0; 
	float refractanceIndex = 1; 
	float absorption = 4;

	bool IsGlass() const
	{
		return type == MaterialType::Dielectric || type == MaterialType::Beers;
	}
};