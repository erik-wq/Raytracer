#pragma once

#include "include/Shapes/Sphere.h"

const int RayChecks = 1000;
const float Speed = 0.02f;

class VoxelGrid;

class Player
{
public:
	Player(Scene* scene);
	~Player();

	void Tick(float delta);

	float3 GetPosition() const
	{
		return sphere->center;
	}

	bool HasMoved() const
	{
		return moved;
	}

	bool HasWon() const
	{
		return won;
	}

private:
	bool CheckMove(const float3& oldPos) const;

	Sphere* sphere;
	Scene* scene;
	VoxelGrid* grid;

	bool moved = false;

	bool won = false;

	float3 velocity = float3(0, -Speed, 0);

	float3* rayDirection;
};