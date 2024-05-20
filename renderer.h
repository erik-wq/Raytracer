#pragma once
#include <vector>
#include "include/Material.h"
#include "include/Lights/DirectionalLight.h"
#include "include/Lights/PointLight.h"
#include "include/Lights/SpotLight.h"

class Skydome;
struct DirectionalLight;
struct PointLight;
struct SpotLight;
class Player;
class Level;

#define DEBUG_INFO 0

namespace Tmpl8
{
	const int MaxBounces = 10;
	const uint8_t maxRayPerPixel = 2;
	const int MaxMaterail = 9;
	const float Epsilon = 0.0001f;
	const float HistoryFilter = 0.9f;
	const float ReprojectionFilter = 0.75f; 

	const float DepthFilter = 0.05f;

class Renderer : public TheApp
{
public:
	// game flow methods
	void Init();
	float3 Trace( Ray& ray, int bounce );
	void Tick( float deltaTime );
	void UI();
	void Shutdown();
	// input handling
	void MouseUp( int  ) { /* implement if you want to detect mouse button presses */ }
	void MouseDown( int  ) { /* implement if you want to detect mouse button presses */ }
	void MouseMove( int x, int y ) { mousePos.x = x, mousePos.y = y; }
	void MouseWheel( float  ) { /* implement if you want to handle the mouse wheel */ }
	void KeyUp( int  ) { /* implement if you want to handle keys */ }
	void KeyDown( int  ) { /* implement if you want to handle keys */ }

	float3 OfsetedRay(const float3& point, const float3& normal) const;

	void ClearLights();

	void StartLevel();

	void BackToMenu();

	SpotLight* AddSpotLight(const SpotLight& light);
	PointLight* AddPointLight(const PointLight& light);
	DirectionalLight* AddDirectionalLight(const DirectionalLight& light);

	// data members
	int2 mousePos;
	float4* accumulator;
	// copy of accumulator for history
	float4* reprojection;

	float* depthBuffer;
	float* reprojectingDepthBuffer;

	Scene scene;
	Camera camera;
private:
	Skydome* sky;
	std::vector<DirectionalLight> directionalLights;
	std::vector<PointLight> pointLights;
	std::vector<SpotLight>	spotLights;
	std::vector<Materail>	materials;

	Level* level;

	uint frameCount = 0;

	int maxBounces = 4;

	bool normals = false;

	float3 SampleLight(const Ray& ray) const;

	float3 CalculateDirectLight(const Ray& ray, const DirectionalLight& Light) const;

	float3 CalculatePointLight(const Ray& ray, const PointLight& light) const;
	
	float3 CalculateSpotLight(const Ray& ray, const SpotLight& light) const;

	Ray ScatterRay(Ray& ray, const float3& normal, const float refractIndex, float3& colorMult, int objectId, bool beer = false, const float absorptionIndex = 0);
	void ResetAccumulator();
	void HandlePointLightUI();
	void HandleSpotLightLightUI();

	void HandleMaterialUI();

	float3 RandomPointSphere(const float3& normal) const;

	float3 VectorAsColor(const float3& vector) const;

	bool Refract(float3& refracted, const float3& normal, const float3& direction, const float n1, const float n2);

	float Fresnell(const Ray& ray, const float3& normal, const float n1, const float n2) const;

	bool Reproject(const Ray& ray, float4& color);

	float3 CheckPosition(const float3& p);
};

}