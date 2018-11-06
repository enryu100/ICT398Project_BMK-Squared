#pragma once
#ifndef _H_MASS_RIGIDBODY_
#define _H_MASS_RIGIDBODY_

#include "Rigidbody.h"
#define GRAVITY_CONST vec3(0.0f, -9.82f, 0.0f)

CollisionManifold FindCollisionFeatures(
		const RigidbodyVolume& ra, const RigidbodyVolume& rb);
void ApplyImpulse(RigidbodyVolume &A,
		RigidbodyVolume &B, const CollisionManifold &M, int c);

class RigidbodyVolume : public Rigidbody {
public:
	vec3 position;
	vec3 velocity;
	vec3 forces; // Sum of all forces
	float mass;
	float cor; //Coefficient of restitution
	float friction;

	vec3 orientation;
	vec3 angVel;
	vec3 torques;//Sum torques

	OBB box;
	Sphere sphere;

	public:
	inline RigidbodyVolume()
		: cor(0.5f), mass(1.0f),
		friction(0.6f) {
				type = RIGIDBODY_TYPE_BASE;
	}
	inline RigidbodyVolume(int bodyType):
		cor(0.5f), mass(1.0f), friction(0.6f){
				type = bodyType;
	}

	~RigidbodyVolume() { }

	void Render();
	void Update(float dt); // UpdatePosition
	void ApplyForces();

	void SynchCollisionVolumes();
	float InvMass();
	void AddLinearImpulse(const vec3 &impulse);
	mat4 InvTensor();
	virtual void AddRotationImpulse(const vec3 &point,
			const vec3 &impulse);

};

#endif
