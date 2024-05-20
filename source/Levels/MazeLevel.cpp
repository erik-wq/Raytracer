#include "precomp.h"
#include "include/Levels/MazeLevel.h"
#include "include/VoxelGrid.h"
#include "include/Player.h"
#include "include/Lights/SpotLight.h"
#include "include/Lights/PointLight.h"

MazeLevel::MazeLevel(Scene* scene, Renderer* renderer) : Level(scene, renderer)
{
	scene->AddVoxelGrid(VoxelGrid());
	player = new Player(scene);

	float3 playerPos = player->GetPosition();
	playerPos.y += 0.4f;

	SpotLight light = { float4(0), playerPos, float3(0,1,0), float3{100 / 255.f,100 / 255.f,100 / 255.f}, 0.95 };
	playerLight = renderer->AddSpotLight(light);

	PointLight finishLight;
	finishLight.position = float3(0.03f, 0.025f, 0.98f);
	finishLight.color = float3(0.05f, 0.05f, 0.05f);
	renderer->AddPointLight(finishLight);

	PointLight firstLigth;
	firstLigth.position = float3(0.79f, 0.025f, 0.6f);
	firstLigth.color = float3(0.05f, 0.05f, 0.05f);
	renderer->AddPointLight(firstLigth);

	UpdateCamera();
}

MazeLevel::~MazeLevel()
{
	delete player;
}

void MazeLevel::Tick(float deltaTime)
{
	if (player->HasWon())
	{
		return;
	}
	player->Tick(deltaTime);
	playerLight->position = player->GetPosition() + float3(0, 0.4f, 0);
	
	if (player->HasMoved())
	{
		UpdateCamera();
	}

	fps = 1.f / deltaTime;
	elapsedTime += deltaTime;
}

void MazeLevel::RenderUI()
{

	if (!player->HasWon())
	{

		ImGui::SetNextWindowPos(ImVec2(0, 0));
		ImGui::SetNextWindowSize(ImVec2(SCRWIDTH / 10, SCRHEIGHT / 10), ImGuiCond_Always);

		ImGui::Begin("New Window", nullptr, ImGuiWindowFlags_NoTitleBar);

		ImGui::Text("Maze Level", ImVec2(SCRWIDTH / 10, SCRHEIGHT / 10));

		ImGui::Text("fps : %.2f", fps , ImVec2(SCRWIDTH / 10, SCRHEIGHT / 10));

		ImGui::Text("time : %.2f", elapsedTime ,ImVec2(SCRWIDTH / 10, SCRHEIGHT / 10));
		
		ImGui::End();
		return;
	}

	ImGui::Begin("New Window", nullptr, ImGuiWindowFlags_NoTitleBar);
	ImGui::SetNextWindowPos(ImVec2(HALFSCRWIDTH - HALFWIDTHWINDOW, HALFSCRHEIGHT - HALFHEIGHTWINDOW), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(HALFWIDTHWINDOW * 2, HALFHEIGHTWINDOW * 2), ImGuiCond_Always);

	ImGui::Begin("New Window", nullptr, ImGuiWindowFlags_NoTitleBar);

	ImGui::SetCursorPos(ImVec2(HALFWIDTHWINDOW / 2, HALFHEIGHTWINDOW / 2));

	if (ImGui::Button("Menu", ImVec2(HALFWIDTHWINDOW, HALFHEIGHTWINDOW / 4)))
	{
		renderer->BackToMenu();
	}

	ImGui::End();
}

void MazeLevel::UpdateCamera() const
{
	float3 playerPos = player->GetPosition();
	float3 cameraPos = playerPos + cameraOffset;
	renderer->camera.UpdateCamera(cameraPos, playerPos);
}
