#pragma once

class Level
{
public:
	Level(Scene* scene, Renderer* renderer) : scene(scene), renderer(renderer) {};
	virtual ~Level() 
	{
		scene->ClearObjects();
		renderer->ClearLights();
	};

	virtual void Tick(float) {};
	virtual void RenderUI() {};

protected:
	Scene* scene;
	Renderer* renderer;
};