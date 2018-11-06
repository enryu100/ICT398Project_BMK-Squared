#pragma once
#ifndef _H_PHYSICS_SYSTEM_
#define _H_PHYSICS_SYSTEM_

#include "Rigidbody.h"

class PhysicsSystem {
protected:
	std ::vector<Rigidbody *> bodies;
	std::vector<OBB> constraints;

	std::vector<Rigidbody *> colliders1;
	std::vector<Rigidbody *> colliders2;
	std::vector<CollisionManifold> results;

	//Try to keep between 0.2 - 0.8
	float LinearProjectionPercent;
	//keep between 0.01 and 0.1
	float PenetrationSlack;
	//[1 to 20], Larger = more accurate
	int ImpulseIteration;

public:
	PhysicsSystem();
	void Update(float deltaTime);
	void Render();
	void AddRigidbody(Rigidbody *body);
	void AddConstraint(const OBB &constratint);

	void ClearRigidbodys();
	void ClearConstraints();
};

#endif
