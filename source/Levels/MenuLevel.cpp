#include "precomp.h"
#include "include/Levels/MenuLevel.h"
#include "include/Shapes/Sphere.h"
#include "include/Shapes/Plane.h"
#include "include/Lights/PointLight.h"

MenuLevel::MenuLevel(Scene* scene, Renderer* renderer) : Level(scene, renderer)
{
	Plane* plane = scene->AddPlane(Plane(float3(0, 0, 0), float3(0, 1, 0)));
	plane->material = 6;

	Sphere* sphere = scene->AddSphere(Sphere(float3(-2, 0.75f, 0), 0.5f));
	sphere->material = 7;
	
	sphere = scene->AddSphere(Sphere(float3(2, 0.75f, 0), 0.5f));
	sphere->material = 8;

	sphere = scene->AddSphere(Sphere(float3(0, 0.75f, 0), 0.5f));
	sphere->material = 9;

	DirectionalLight light;
	light.direction = normalize(float3(-0.7f, -0.75f,0));
	light.color = float3(1, 1, 1);
	renderer->AddDirectionalLight(light);

	renderer->camera.UpdateCamera(float3(0, 1.2f, -4), float3(0, 1, 0));
}

void MenuLevel::RenderUI()
{
	if (!cameraSet)
	{
		renderer->camera.UpdateCamera(float3(0, 1.2f, -4), float3(0, 1, 0));
		cameraSet = true;
	}

	ImGui::SetNextWindowPos(ImVec2(HALFSCRWIDTH - HALFWIDTHWINDOW, HALFSCRHEIGHT - (HALFHEIGHTWINDOW / 2)), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(HALFWIDTHWINDOW * 2, HALFHEIGHTWINDOW), ImGuiCond_Always);

	ImGui::Begin("New Window", nullptr, ImGuiWindowFlags_NoTitleBar);

	ImGui::SetCursorPos(ImVec2(HALFWIDTHWINDOW / 2, HALFHEIGHTWINDOW / 4));
	if (ImGui::Button("Play", ImVec2(HALFWIDTHWINDOW, HALFHEIGHTWINDOW / 4)))
	{
		renderer->StartLevel();
	}

	ImGui::End();
}
