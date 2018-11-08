#pragma once
#ifndef _H_RIGIDBODY_
#define _H_RIGIDBODY_

#define RIGIDBODY_TYPE_BASE 0
#define RIGIDBODY_TYPE_PARTICLE 1
#define RIGIDBODY_TYPE_SPHERE 2
#define RIGIDBODY_TYPE_BOX 3

#include "Geometry3d.h"
#include <vector>

class Rigidbody {
public:
	int type;
	Rigidbody();
	virtual ~Rigidbody() {}

	virtual void Update(float deltaTime) {}
	virtual void Render() {}
	virtual void ApplyForces() {}
	virtual void SolveConstratins(
			const std::vector<OBB> &constraints) {}
	inline Rigidbody() {
		type = RIGIDBODY_TYPE_BASE;
	}

	inline bool HasVolume() {
		return type == RIGIDBODY_TYPE_SPHERE ||
			   type == RIGIDBODY_TYPE_BOX;
	}
};

#endif
