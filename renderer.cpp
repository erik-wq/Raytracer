#include "precomp.h"
#include "include/Skydome.h"
#include "include/Shapes/Sphere.h"
#include "include/VoxelGrid.h"
#include "include/Shapes/Plane.h"
#include "include/Levels/MazeLevel.h"
#include "include/Levels/MenuLevel.h"

// REFERENCE IMAGES:
// https://www.rockpapershotgun.com/minecraft-ray-tracing
// https://assetsio.reedpopcdn.com/javaw_2019_04_20_23_52_16_879.png
// https://www.pcworld.com/wp-content/uploads/2023/04/618525e8fa47b149230.56951356-imagination-island-1-on-100838323-orig.jpg

// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
	InitSeed(static_cast<uint>(time(nullptr)));

	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);

	reprojection = (float4*)MALLOC64(SCRWIDTH * SCRHEIGHT * 16);
	memset(reprojection, 0, SCRWIDTH * SCRHEIGHT * 16);

	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);

	depthBuffer = (float*)MALLOC64(SCRWIDTH * SCRHEIGHT * 4);
	memset(depthBuffer, 0, SCRWIDTH * SCRHEIGHT * 4);

	reprojectingDepthBuffer = (float*)MALLOC64(SCRWIDTH * SCRHEIGHT * 4);
	memset(reprojectingDepthBuffer, 0, SCRWIDTH * SCRHEIGHT * 4);
	
	// try to load a camera
	FILE* f = fopen("camera.bin", "rb");
	if (f)
	{
		fread(&camera, 1, sizeof(Camera), f);
		fclose(f);
	}

	level = new MenuLevel(&scene, this);
	sky = new Skydome();
	materials = std::vector<Materail>(MaxMaterail);
	materials[3].albedo = {0,0,1};
	materials[2].albedo = {1,0,0};
	materials[1].albedo = {0,1,0};
	materials[0].albedo = {1,0,0};

	materials[4].albedo = {0,0,1};

	materials[5].albedo = {0.5f,0.5f,0.5f};
	materials[5].reflectance = 0.25f;

	materials[6].albedo = {1,0,0};
	materials[6].reflectance = 1;

	materials[7].type = MaterialType::Beers;
	materials[7].refractanceIndex = 1.f;
	materials[7].absorption = 1;
	materials[7].AbsorbedColor = float3(0.9f, 0, 0);

	materials[8].type = MaterialType::Dielectric;
	materials[8].refractanceIndex = 1.56f;
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Trace(Ray& ray, int bounce)
{

	int ID = scene.FindNearest(ray);

	if (ray.voxel == 0)
	{
		return  sky->GetSkyColor(ray);
	}

	if (bounce >= MaxBounces)
	{
		// printf("reaching rec limit \n");
		return { 0 };
	}

	float3 I = ray.IntersectionPoint();

	float3 normal = ray.GetNormal();

	if (normals)
	{
		return VectorAsColor(normal);
	}

	Materail material = materials[ray.voxel - 1];

	if (material.IsGlass())
	{
		bool beer = material.type == MaterialType::Beers;
		float3 colorMult = float3(1);
		Ray newRay = ScatterRay(ray, normal, material.refractanceIndex, colorMult, ID, beer, material.absorption);
		float3 color = Trace(newRay, bounce + 1);
		if (beer)
		{
			float3 abs = (1 - colorMult) * material.AbsorbedColor;
			return color * colorMult + abs;
		}
		return color;
	}

	float3 reflection = float3(0);

	if (material.reflectance > 0 && bounce < maxBounces)
	{
		// calculate reflection ray
		float3 reflectedDirection = normalize(reflect(normalize(ray.D), normal));

		Ray reflectionRay(OfsetedRay(I, normal), reflectedDirection);

		if (material.reflectance == 1)
		{
			return Trace(reflectionRay, bounce + 1);
		}

		reflection = Trace(reflectionRay, bounce + 1);
	}

	float3 albedo = material.albedo * SampleLight(ray);
	return (albedo * (1 - material.reflectance)) + reflection * material.reflectance;
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{

	// pixel loop
	Timer t;

	memcpy(reprojection, accumulator, SCRWIDTH * SCRHEIGHT * 16);
	memcpy(reprojectingDepthBuffer, depthBuffer, SCRWIDTH * SCRHEIGHT * 4);

	float samplesMult = 1.0f / maxRayPerPixel;

	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++)
		{
			const int pos = x + y * SCRWIDTH;

			Ray depthRay = camera.GetPrimaryRay(static_cast<float>(x), static_cast<float>(y));
			Trace(depthRay, MaxBounces - 1);

			float3 position = depthRay.IntersectionPoint();
			
			float depthdiff = depthRay.t - depthBuffer[pos];

			const bool changeDepth = !isnan(depthdiff) && depthdiff != 0;

			depthBuffer[pos] = depthRay.t;

			float3 totalLight{ 0 };

			for (int i = 0; i < maxRayPerPixel; i++)
			{
				float randX = Rand(1);
				float randY = Rand(1);
				Ray primaryRay = camera.GetPrimaryRay(static_cast<float>(x) + randX, static_cast<float>(y) + randY);

				totalLight += Trace(primaryRay, 0);
			}

			const float4 newPixel = float4{ totalLight * samplesMult, 0.0f };

			if (camera.moved)
			{
				float4 col;

				if (Reproject(depthRay, col))
				{
					accumulator[pos] = col * ReprojectionFilter + (1 - ReprojectionFilter) * newPixel;
				}
				else
				{
					accumulator[pos] = newPixel;
				}
			}
			else
			{

				if (changeDepth)
				{
					accumulator[pos] = newPixel;
				}
				else
				{
					accumulator[pos] = accumulator[pos] = accumulator[pos] * HistoryFilter + (1 - HistoryFilter) * newPixel;
				}
			}
			float4 Color = accumulator[pos];
			screen->pixels[pos] = RGBF32_to_RGB8(&Color);
		}
	}

	if (camera.moved)
	{
		camera.moved = false;
	}

#if DEBUG_INFO
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	printf("%5.2fs \n", avg / 1000);
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000.0f / avg, rps = (SCRWIDTH * SCRHEIGHT) / avg;
	printf("%5.2fms (%.1ffps) - %.1fMrays/s\n", avg, fps, rps / 1000);
	// handle user input
#endif
	camera.HandleInput(deltaTime);
	level->Tick(deltaTime * 0.001f);
}


void Tmpl8::Renderer::ResetAccumulator()
{
	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
	frameCount = 1;
}

float3 Tmpl8::Renderer::VectorAsColor(const float3& vector) const
{
	return normalize(vector) * 0.5f + 0.5f;
}

#pragma region LightCalculations

float3 Tmpl8::Renderer::SampleLight(const Ray& ray) const
{
	const int point = (int)pointLights.size();
	const int spot = (int)spotLights.size();
	const int direct = (int)directionalLights.size();
	const int sum = point + spot + direct;

	int x = (int)(ceilf(Rand((float)sum))) - 1;

	if (x < 0)
	{
		return { 0 };
	}

	// only one direct per level
 	if (direct > 0 && x < direct)
	{
		return CalculateDirectLight(ray, directionalLights[x]);
	}
	x -= direct;

	if (point > 0 && x < point)
	{
		return CalculatePointLight(ray, pointLights[x]);
	}
	x -= point;

	if (spot > 0)
	{
		return CalculateSpotLight(ray, spotLights[x]);
	}

	return 0;
}
float3 Tmpl8::Renderer::CalculateDirectLight(const Ray& ray, const DirectionalLight& Light) const
{
	float intes = Light.color.x + Light.color.y + Light.color.z;
	if (intes <= 0.001)
	{
		return 0;
	}

	const float3 intersectionPoint = ray.IntersectionPoint();
	const float3 dir = -Light.direction;

	const float3 normal = ray.GetNormal();
	//light angle
	const float cosTheta = dot(dir, normal);
	if (cosTheta <= 0)
		return 0;

	const float3 lightIntensity = max(0.0f, cosTheta) * Light.color;
	Ray shadowRay(OfsetedRay(intersectionPoint, normal), dir);

	if (scene.IsOccluded(shadowRay))
		return 0;

	return lightIntensity;
}

float3 Tmpl8::Renderer::CalculatePointLight(const Ray& ray, const PointLight& light) const
{
	float intes = light.color.x + light.color.y + light.color.z;
	if (intes <= 0.001)
	{
		return 0;
	}

	const float3 intersectionPoint = ray.IntersectionPoint();
	const float3 dir = light.position - intersectionPoint;
	const float dst = length(dir);
	const float3 dirNormalized = dir / dst;

	const float3 normal = ray.GetNormal();
	//light angle
	const float cosTheta = dot(dirNormalized, normal);
	if (cosTheta <= 0)
		return 0;
	const float3 lightIntensity = max(0.0f, cosTheta) * light.color / (dst * dst);

	auto point = OfsetedRay(light.position, normal);
	point = OfsetedRay(point, normal);

	Ray shadowRay(point, -dirNormalized);
	shadowRay.t = dst;
	if (scene.IsOccluded(shadowRay))
		return 0;

	return lightIntensity;
}

float3 Tmpl8::Renderer::CalculateSpotLight(const Ray& ray, const SpotLight& light) const
{
	float intes = light.color.x + light.color.y + light.color.z;
	if (intes <= 0.001)
	{
		return 0;
	}

	const float3 intersectionPoint = ray.IntersectionPoint();
	const float3 dir = light.position - intersectionPoint;
	float dst = length(dir);
	float3 dirNormalized = dir / dst;

	float3 normal = ray.GetNormal();
	//light angle
	float cosTheta = dot(dirNormalized, light.direction);
	if (cosTheta <= light.angle)
		return 0;

	float alphaCutOff = 1.0f - (1.0f - cosTheta) * 1.0f / (1.0f - light.angle);

	const float distRcp = 1 / (dst * dst);
	float3 lightIntensity = max(0.0f, cosTheta) * light.color * distRcp;

	Ray shadowRay(OfsetedRay(light.position, normal), -dirNormalized);
	shadowRay.t = dst;
	if (scene.IsOccluded(shadowRay))
		return 0;

	return lightIntensity * alphaCutOff;
}
#pragma endregion

#pragma region RayCalculations
float3 Tmpl8::Renderer::OfsetedRay(const float3& point, const float3& normal) const
{
	return point + normal * Epsilon;
}

void Tmpl8::Renderer::ClearLights()
{
	pointLights.clear();
	spotLights.clear();
	directionalLights.clear();
}

void Tmpl8::Renderer::StartLevel()
{
	delete level;
	level = new MazeLevel(&scene, this);
	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
}

void Tmpl8::Renderer::BackToMenu()
{
	delete level;
	level = new MenuLevel(&scene, this);
	memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
}

SpotLight* Tmpl8::Renderer::AddSpotLight(const SpotLight& light)
{
	spotLights.push_back(light);
	return &spotLights.back();
}

PointLight* Tmpl8::Renderer::AddPointLight(const PointLight& light)
{
	pointLights.push_back(light);
	return &pointLights.back();
}

DirectionalLight* Tmpl8::Renderer::AddDirectionalLight(const DirectionalLight& light)
{
	directionalLights.push_back(light);
	return &directionalLights.back();
}

Ray Tmpl8::Renderer::ScatterRay(Ray& ray, const float3& normal, const float refractIndex, float3& colorMult, int objectId, bool beer,  const float absorptionIndex)
{
	float3 dir;
	float3 I = ray.IntersectionPoint();
	if (Refract(dir, normal, ray.D, 1, refractIndex))
	{
		Ray refractedRay(I + dir * 0.001f, dir);
		refractedRay.type = ray.type;

		scene.FindExit(refractedRay, objectId);

		if (beer)
		{
			const float absorption = (absorptionIndex / 2 )* refractedRay.t;
			const float mult = exp(-absorption);
			// printf("dist %lf, mult %lf\n",refractedRay.t, mult);
			colorMult = float3(mult, mult, mult);
		}

		I = refractedRay.IntersectionPoint();
		Refract(dir, -refractedRay.GetNormal(), refractedRay.D, refractIndex, 1);
		return Ray(I + dir * 0.0001f, dir);
	}

	dir = reflect(ray.D, normal);
	return Ray(I + dir * 0.0001f, dir);
}

float3 Tmpl8::Renderer::RandomPointSphere(const float3& normal) const
{
	while (true)
	{
		// Generate a random point within a cube
		const float3 point{ RandomFloat(), RandomFloat(), RandomFloat() };
		// If the point is within the unit sphere
		if (sqrLength(point) < 1)
		{
			// Normalize the point
			// Use it as a unit vector
			const float3 vec{ normalize(point) };

			// Return Lambertian reflection vector
			return normal + vec;
		}
	}
}

bool Tmpl8::Renderer::Refract(float3& refracted, const float3& normal, const float3& direction, const float n1, const float n2)
{
	float eta = n1 / n2;
	float cos = dot(normal, -direction);
	float w = eta * cos;

	float k = 1.f + (w - eta) * (w + eta);
	if (k < 0)
	{
		return false;
	}

	float3 t = (w - sqrt(k)) * normal - eta * (-direction);

	refracted = normalize(t);
	return true;
}

float Tmpl8::Renderer::Fresnell(const Ray& ray, const float3& normal, const float n1, const float n2) const
{
	const float eta = n1 / n2;

	const float cosAi = dot(normal, -ray.D);

	const float s = eta * sin(cosAi);

	const float cosAt = sqrtf(1 - (s * s));

	const float n1CosAi = n1 * cosAi;
	const float n2CosAi = n2 * cosAi;
	const float n1CosAt = n1 * cosAt;
	const float n2CosAt = n2 * cosAt;

	float a = (n1CosAi - n2CosAt) * 1 / (n1CosAi + n2CosAt);
	float b = (n1CosAt - n2CosAi) * 1 / (n1CosAt + n2CosAi);

	return 0.5f * (a * a + b * b);
}

#pragma endregion

#pragma region UI

// -----------------------------------------------------------
// Update user interface (imgui)
// -----------------------------------------------------------
void Renderer::UI()
{
#if DEBUG_INFO
	// ray query on mouse
	ImGui::DragFloat3("camera position : ", camera.camPos.cell);
	ImGui::DragFloat3("camera target : ", camera.camTarget.cell);

	Ray r = camera.GetPrimaryRay((float)mousePos.x, (float)mousePos.y);
	scene.FindNearest(r);
	ImGui::Begin("Info");
	ImGui::Text("voxel: %i", r.voxel);

	ImGui::SliderInt("reflection bounces :", &maxBounces, 0, MaxBounces);

	HandleMaterialUI();

	/*
	if (ImGui::CollapsingHeader("Direct Lights"))
	{
		ImGui::ColorEdit3("Direct Color:", directional.color.cell);
		ImGui::SliderFloat3("Direct Direction:", directional.direction.cell, -1.0f, 1.0f);
		directional.direction = normalize(directional.direction);
	}
	*/

	HandlePointLightUI();
	HandleSpotLightLightUI();

	ImGui::Checkbox("normals:", &normals);

	ImGui::End();
#endif
	level->RenderUI();
}

void Tmpl8::Renderer::HandlePointLightUI()
{
	if (!ImGui::CollapsingHeader("Point Lights"))
	{
		return;
	}

	for (int i = 0; i < pointLights.size(); i++)
	{
		if (!ImGui::CollapsingHeader(("Point Lights " + to_string(i)).c_str()))
		{
			continue;
		}
		ImGui::ColorEdit3(("Point Color: " + to_string(i)).c_str(), pointLights[i].color.cell);
		ImGui::SliderFloat3(("Point Position: " + to_string(i)).c_str(), pointLights[i].position.cell, 0.f, 2.f);
	}
}

void Tmpl8::Renderer::HandleSpotLightLightUI()
{
	if (!ImGui::CollapsingHeader("SpotLight Lights"))
	{
		return;
	}

	for (int i = 0; i < spotLights.size(); i++)
	{
		if (!ImGui::CollapsingHeader(("SpotLight Lights " + to_string(i)).c_str()))
		{
			continue;
		}
		ImGui::ColorEdit3(("SpotLight Color: " + to_string(i)).c_str(), spotLights[i].color.cell);
		ImGui::SliderFloat3(("SpotLight Position: " + to_string(i)).c_str(), spotLights[i].position.cell, -10.0f, 10.0f);
		ImGui::SliderFloat3(("SpotLight Direction: " + to_string(i)).c_str(), spotLights[i].direction.cell, -1.0f, 1.0f);
		ImGui::SliderFloat(("SpotLight Angle: " + to_string(i)).c_str(), &spotLights[i].angle, -1.0f, 1.0f);
		spotLights[i].direction = normalize(spotLights[i].direction);
	}
}
void Tmpl8::Renderer::HandleMaterialUI()
{
	if (ImGui::CollapsingHeader("Material"))
	{
		int i = 0;
		for (auto& Mat : materials)
		{
			ImGui::ColorEdit3(("Material Color:" + to_string(i)).c_str(), Mat.albedo.cell);
			ImGui::SliderFloat(("Material Reflection:" + to_string(i)).c_str(), &Mat.reflectance, 0, 1.0f);
			ImGui::SliderFloat(("Material Refractance:" + to_string(i)).c_str(), &Mat.refractanceIndex, 1, 2.5f);
			ImGui::SliderFloat(("Material Absorption:" + to_string(i)).c_str(), &Mat.absorption, 0.f, 3.0f);
			ImGui::ColorEdit3(("Material Absorption:" + to_string(i)).c_str(), Mat.AbsorbedColor.cell);
			i++;
		}
	}
}
#pragma endregion


// -----------------------------------------------------------
// User wants to close down
// -----------------------------------------------------------
void Renderer::Shutdown()
{
	// save current camera
	FILE* f = fopen("camera.bin", "wb");
	fwrite(&camera, 1, sizeof(Camera), f);
	fclose(f);
	delete sky;
	delete level;
}

// add depth rejection
bool Tmpl8::Renderer::Reproject(const Ray& ray, float4& color)
{
	// material 4 is player he is moving all the ti 
	if (ray.t == 1e34f || ray.voxel == 4)
	{
		return false;
	}

	const float3 point = ray.IntersectionPoint();
	float3 reprojectedUV = CheckPosition(point);
	float2 uv = float2{ reprojectedUV.x, reprojectedUV.y };
	uv.x /= camera.aspect;
	uv = uv * 0.5f + 0.5f;
	uv.y = 1 - uv.y;

	if (uv.x > 0.0 && uv.x < 1.0 && uv.y > 0.0 && uv.y < 1.0 && reprojectedUV.z > -0.5f)
	{
		uv.x *= SCRWIDTH;
		uv.y *= SCRHEIGHT;

		int pos = int(floor(uv.x) + floor(uv.y) * SCRWIDTH);
		float oldDepth = reprojectingDepthBuffer[pos];

		if (oldDepth == 1e34f)
		{
			return false;
		}

		float3 oldDirection = camera.GetOldDirection(uv.x, uv.y);
		float3 oldPoint = camera.oldPosition + oldDirection * oldDepth;

		float dist = length(oldPoint - point);
		if (dist < 0.009f)
		{
			color = reprojection[pos];
			return true;
		}
	}
	return false;
}

float3 Tmpl8::Renderer::CheckPosition(const float3& p)
{
	float3 toPoint = p - camera.oldPosition;

	float3 toPointNrm = normalize(toPoint);
	float3 direction = normalize(camera.oldTarget - camera.oldPosition);

	float3 right = normalize(cross(float3(0.0, 1.0, 0), direction));
	float3 up = (cross(direction, right));

	float d = dot(direction, toPointNrm);
	if (d < 0.01)
	{
		return float3(0.0, 0.0, -1.0);
	}

	float3 fwd = direction * 2;
	d = 2 / d;

	toPointNrm = toPointNrm * d - fwd;

	float x = dot(toPointNrm, right);
	float y = dot(toPointNrm, up);

	return float3{ x, y, 1.0 };
}
