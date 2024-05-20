#pragma once

class Light
{
public:
	virtual ~Light() = default;
	/**
	* multiplay with material color
	*/
	virtual float3 CalculateLight(Ray& ray, Scene& scene, Renderer* renderer);
};