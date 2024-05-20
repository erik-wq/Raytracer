#pragma once
#include "Level.h"

class Player;

class MazeLevel : public Level
{
public:
	MazeLevel(Scene* scene, Renderer* renderer);
	~MazeLevel();

	virtual void Tick(float deltaTime) override;
	virtual void RenderUI() override;
private:
	void UpdateCamera() const;

	float3 cameraOffset = float3(0, 0.35f, -0.05f);

	SpotLight* playerLight;
	Player* player;

	float fps;
	float elapsedTime = 0;
};