#include "FixedFunctionPrimitives.h"
#include "PhysicsSystem.h"
#include "RigidbodyVolume.h"
#include "glad/glad.h"

PhysicsSystem::PhysicsSystem() {
	LinearProjectionPercent = 0.45f;
	PenetrationSlack = 0.01f;
	ImpulseIteration = 5;

	colliders1.reserve(100);
	colliders2.reserve(100);
	results.reserve(100);
}

void PhysicsSystem::AddRigidbody(Rigidbody *body) {
	bodies.push_back(body);
}

void PhysicsSystem::AddConstraint(const OBB &obb) {
	constraints.push_back(obb);
}

void PhysicsSystem::ClearRigidbodys() {
	bodies.clear();
}

void PhysicsSystem::ClearConstraints() {
	constraints.clear();
}

void PhysicsSystem::Render() {
	static const float rigidbodyDiffuse[]{
		200.0f / 255.0f, 0.0f, 0.0f, 0.0f
	};
	static const float rigidbodyAmbient[]{
		200.0f / 255.0f, 50.0f / 255.0f, 50.0f / 255.0f, 0.0f
	};
	static const float constraintDiffuse[]{
		0.0f, 200.0f / 255.0f, 0.0f, 0.0f
	};
	static const float constraintAmbient[]{
		50.0f / 255.0f, 200.0f / 255.0f, 50.0f / 255.0f, 0.0f
	};
	static const float zero[] = { 0.0f, 0.0f, 0.0f, 0.0f };

	glColor3f(rigidbodyDiffuse[0],
			rigidbodyDiffuse[1],
			rigidbodyDiffuse[2]);
	glLightfv(GL_LIGHT0, GL_AMBIENT, rigidbodyAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, rigidbodyDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, zero);

	for (int i = 0, size = bodies.size(); i < size; ++i) {
		bodies[i]->Render();
	}

	glColor3f(constraintDiffuse[0],
			constraintDiffuse[1],
			constraintDiffuse[2]);
	glLightfv(GL_LIGHT0, GL_AMBIENT, constraintAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, constraintDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, zero);

	for (int i = 0; i < constraints.size(); ++i) {
		::Render(constraints[i]);
	}
}

void PhysicsSystem::Update(float deltaTime) {

	colliders1.clear();
	colliders2.clear();
	results.clear();

	for (int i = 0, size = bodies.size(); i < size; ++i) {
		bodies[i]->ApplyForces();
		for (int j = i; j < size; ++j) {
			if (i == j) {
				continue;
			}
			CollisionManifold result;
			ResetCollisionManifold(&result);
			if (bodies[i]->HasVolume() &&
					bodies[j]->HasVolume()) {
				RigidbodyVolume *m1 =
						(RigidbodyVolume *)bodies[i];
				RigidbodyVolume *m2 =
						(RigidbodyVolume *)bodies[j];
				result = FindCollisionFeatures(*m1, *m2);
			}
			if (result.colliding) {
				colliders1.push_back(bodies[i]);
				colliders2.push_back(bodies[j]);
				results.push_back(result);
			}
		}
	}
	//Calculate forces acting on the object
	for (int i = 0, size = bodies.size(); i < size; ++i) {
		bodies[i]->ApplyForces();
	}
	for (int k = 0; k < ImpulseIteration; ++k) {
		for (int i = 0; i < results.size(); ++i) {
			int jSize = results[i].contacts.size();
			for (int j = 0; j < jSize; ++j) {
				if (colliders1[i]->HasVolume() && colliders2[i]->HasVolume()) {
					RigidbodyVolume* m1 =
							(RigidbodyVolume*)colliders1[i];
					RigidbodyVolume* m2 =
							(RigidbodyVolume*)colliders2[i];
					ApplyImpulse(*m1, *m2, results[i], j);
				}
			}
		}
	}
	for (int i = 0, size = bodies.size(); i < size; ++i) {
		bodies[i]->Update(deltaTime);
	}

	for (int i = 0, size = bodies.size(); i < size; ++i) {
		if (!colliders1[i]->HasVolume() && !colliders2[i]->HasVolume()){
			RigidbodyVolume *m1 = (RigidbodyVolume *)colliders1[i];
			RigidbodyVolume *m2 = (RigidbodyVolume *)colliders2[i];
			float totalMass = m1->InvMass() + m2->InvMass();

			if (totalMass == 0.0f) {
				continue;
			}
			float depth = fmaxf(results[i].depth - PenetrationSlack, 0.0f);
			float scalar = depth / totalMass;
			vec3 correction = results[i].normal * scalar * LinearProjectionPercent;

			m1->position = m1->position - correction * m1->InvMass();
			m2->position = m2->position + correction * m2->InvMass();

			m1->SynchCollisionVolumes();
			m2->SynchCollisionVolumes();
		}
	}

	for (int i = 0, size = bodies.size(); i < size; ++i) {
		bodies[i]->SolveConstratins(constraints);
	}
}
