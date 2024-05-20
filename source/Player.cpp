#include "precomp.h"
#include "include/Player.h"
#include "include/VoxelGrid.h"
#include <set>

Player::Player(Scene* scene) : scene(scene)
{
	float voxelSize = 1.f / GRIDSIZE;

	sphere = scene->AddSphere(Sphere(float3(0.015f, 0.02f, 0.015f), voxelSize));
	sphere->material = 4;

	grid = scene->GetVoxelGrid(0);

	rayDirection = new float3[RayChecks];

	float phi = PI * (3.f - sqrt(5.f));

	// generate equaly distributes direction rays for sphere
	for (int i = 0; i < RayChecks; i++)
	{
		float y = 1 - (i / float(RayChecks - 1)) * 2;
		float radius = sqrt(1 - y * y);

		float theta = phi * i;

		float x = cos(theta) * radius;
		float z = sin(theta) * radius;
		float3 dir = normalize(float3(x,y,z));
		rayDirection[i] = dir;
	}
}

Player::~Player()
{
	delete[] rayDirection;
}

void Player::Tick(float delta)
{
	if (won)
	{
		return;
	}

	moved = false;

	if (IsKeyDown(GLFW_KEY_S))
	{
		velocity.z = -Speed;
	}
	else if (IsKeyDown(GLFW_KEY_W))
	{
		velocity.z = Speed;
	}
	else
	{
		velocity.z = 0;
	}

	if (IsKeyDown(GLFW_KEY_A))
	{
		velocity.x = -Speed;
	}
	else if (IsKeyDown(GLFW_KEY_D))
	{
		velocity.x = Speed;
	}
	else
	{
		velocity.x = 0;
	}

	if (velocity.x == 0 && velocity.z == 0)
	{
		return;
	}

	float3 oldPos = sphere->center;

	float3 deltaMove = float3(0, 0, 0);

	float3 move = velocity * delta;
	const float lengt = dot(move, move);
	const int steps = (int)ceilf(lengt / sphere->radius);
	const float multiplier = 1.0f / steps;
	
	int hitCount = 0;

	for (int n = 0; n < steps; n++)
	{
		sphere->center += move * multiplier;
		// check for collision
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < RayChecks; i++)
		{
			Ray ray = Ray(sphere->center, rayDirection[i]);
			ray.t = sphere->radius;

			int3 pos;
			grid->Traverse(ray, pos);

			if (ray.t < sphere->radius)
			{
				if (ray.voxel == 5)
				{
					grid->RemoveVoxel(pos);
					continue;
				}
				velocity.y = 0;
				hitCount++;
				// stop falling falling only at start 
				deltaMove += -ray.D * (sphere->radius - ray.t);
			}
		}

		if (hitCount <= 0)
		{
			break;
		}

		deltaMove *= 1.0f / hitCount;
		sphere->center += deltaMove;
	}


	Ray ray = Ray(sphere->center, float3(0,-1,0));
	ray.t = sphere->radius + 0.005f;

	int3 pos;
	grid->Traverse(ray, pos);

	if (ray.voxel == 3)
	{
		won = true;
	}

	moved = CheckMove(oldPos);
}

bool Player::CheckMove(const float3& oldPos) const
{
	float3 diff = sphere->center - oldPos;
	return diff.x != 0 || diff.y != 0 || diff.z != 0;
}
