#pragma once
#include "Level.h"

class MenuLevel : public Level
{
public:
	MenuLevel(Scene* scene, Renderer* renderer);
	virtual void RenderUI() override;
private:
	bool cameraSet = false;
};