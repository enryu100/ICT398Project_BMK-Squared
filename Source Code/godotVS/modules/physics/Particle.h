#pragma once
#ifndef _H_PARTICLE_
#define _H_PARTICLE_

#include "Rigidbody.h"

class Particle : public Rigidbody {
	vec3 position, oldPosition;
	vec3 forces;
	float mass, bounce;
	vec3 gravity;
	float friction;

public:
	Particle();

	void Update(float deltaTime);
	void Render();
	void ApplyForces();
	void SolveConstratins(
			const std::vector<OBB> &constraints);

	void SetPosition(const vec3 &pos);
	vec3 GetPosition();
	void SetBounce(float b);
	float GetBounce();
};

#endif
