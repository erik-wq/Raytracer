#pragma once

// default screen resolution
#define SCRWIDTH 1024
#define SCRHEIGHT 640

#define HALFSCRWIDTH	SCRWIDTH / 2
#define HALFSCRHEIGHT	SCRHEIGHT / 2

#define	HALFWIDTHWINDOW SCRWIDTH / 8
#define HALFHEIGHTWINDOW SCRHEIGHT / 8


// #define FULLSCREEN
// #define DOUBLESIZE

#define WIDTHMLT 1 / SCRWIDTH
#define HEIGHTMLT 1 / SCRHEIGHT

namespace Tmpl8 {

class Camera
{
public:
	Camera()
	{
		// setup a basic view frustum
		camPos = float3( 0, 0, -2 );
		camTarget = float3( 0, 0, -1 );
		topLeft = float3( -aspect, 1, 0 );
		topRight = float3( aspect, 1, 0 );
		bottomLeft = float3( -aspect, -1, 0 );
		forward = normalize(camTarget - camPos);

		camPosSSE = _mm_set_ps(0, camPos.z, camPos.y, camPos.x);

		widthDiff = topRight - topLeft;
		heightDiff = bottomLeft - topLeft;
	}

	void UpdateCamera(const float3& newPos, const float3& newTarget)
	{


		oldPosition = camPos;
		oldTarget = camTarget;
		oldUp = up;
		oldRight = right;

		oldBottomLeft = bottomLeft;
		oldTopLeft = topLeft;
		oldTopRight = topRight;

		camPos = newPos;
		camTarget = newTarget;

		forward = normalize(camTarget - camPos);
		up = normalize(cross(forward, right));
		right = normalize(cross(up, forward));
		topLeft = camPos + 2 * forward - aspect * right + up;
		topRight = camPos + 2 * forward + aspect * right + up;
		bottomLeft = camPos + 2 * forward - aspect * right - up;

		widthDiff = topRight - topLeft;
		heightDiff = bottomLeft - topLeft;

		camPosSSE = _mm_set_ps(0, camPos.z, camPos.y, camPos.x);

		moved = true;
	}

	float3 GetOldDirection(const float x, const float y)
	{
		// calculate pixel position on virtual screen plane
		const float u = x * (1.0f / SCRWIDTH);
		const float v = y * (1.0f / SCRHEIGHT);
		const float3 P = oldTopLeft + u * (oldTopRight - oldTopLeft) + v * (oldBottomLeft - oldTopLeft);
		// return Ray( oldPosition, normalize( P - oldPosition ) );
		return normalize(P - oldPosition);
		// Note: no need to normalize primary rays in a pure voxel world
	}
	Ray GetPrimaryRay( const float x, const float y )
	{
		// calculate pixel position on virtual screen plane
		const float3 u = x * WIDTHMLT * widthDiff;
		const float3 v = y * HEIGHTMLT * heightDiff;
		const float3 P = topLeft + u + v;

		float3 dir = normalize(P - camPos);

		return Ray( camPos, dir);
	}
	bool HandleInput( const float t )
	{
		if (!WindowHasFocus()) return false;
		oldPosition = camPos;
		oldTarget = camTarget;
		oldUp = up;
		oldRight = right;

		oldBottomLeft = bottomLeft;
		oldTopLeft = topLeft;
		oldTopRight = topRight;

		float speed = 0.0025f * t;
		float3 ahead = normalize( camTarget - camPos );
		float3 tmpUp( 0, 1, 0 );
		right = normalize( cross( tmpUp, ahead ) );
		up = normalize( cross( ahead, right ) );
		bool changed = false;
		if (IsKeyDown( GLFW_KEY_R )) camPos += speed * 2 * up, changed = true;
		if (IsKeyDown( GLFW_KEY_F )) camPos -= speed * 2 * up, changed = true;
		camTarget = camPos + ahead;
		if (IsKeyDown( GLFW_KEY_UP )) camTarget -= speed * up, changed = true;
		if (IsKeyDown( GLFW_KEY_DOWN )) camTarget += speed * up, changed = true;
		if (IsKeyDown( GLFW_KEY_LEFT )) camTarget -= speed * right, changed = true;
		if (IsKeyDown( GLFW_KEY_RIGHT )) camTarget += speed * right, changed = true;
		if (!changed) return false;
		forward = normalize( camTarget - camPos );
		up = normalize( cross( ahead, right ) );
		right = normalize( cross( up, ahead ) );
		topLeft = camPos + 2 * forward - aspect * right + up;
		topRight = camPos + 2 * forward + aspect * right + up;
		bottomLeft = camPos + 2 * forward - aspect * right - up;

		widthDiff = topRight - topLeft;
		heightDiff = bottomLeft - topLeft;

		camPosSSE = _mm_set_ps(0, camPos.z, camPos.y, camPos.x);

		moved = true;
		return true;
	}
	float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
	float3 camPos, camTarget;
	float3 topLeft, topRight, bottomLeft;

	__m128 camPosSSE;

	float3 forward;
	float3 right;
	float3 up;

	float3 oldPosition;
	float3 oldTarget;
	float3 oldUp;
	float3 oldRight;

	float3 widthDiff, heightDiff;

	float3 oldTopLeft, oldTopRight, oldBottomLeft;
	bool moved = false;
};

}