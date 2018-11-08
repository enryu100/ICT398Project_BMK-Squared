#pragma once
#ifndef _H_GEOMETRY_3D_
#define _H_GEOMETRY_3D_

#include "MathVectors.h"
#include "Matrices.h"
#include <vector>

#define AABBSphere(aabb, sphere) \
	SphereAABB(sphere, aabb)

#define OBBSphere(obb, sphere) \
	SphereOBB(sphere, obb)

#define PlaneSphere(plane, sphere) \
	SpherePlane(sphere, plane)

#define OBBAABB(obb, aabb) \
	AABBOBB(aabb, obb)

#define PlaneAABB(plane, aabb) \
	AABBPlane(aabb, plane)

#define PlaneOBB(plane, obb) \
	OBBPlane(obb, plane)

#define SphereTriangle(s, t) \
	TriangleSphere(t, s)

#define AABBTriangle(a, t) \
	TriangleAABB(t, a)

#define OBBTriangle(o, t) \
	TriangleOBB(t, o)

#define PlaneTriangle(p, t) \
	TrianglePlane(t, p)

typedef vec3 Point;

typedef struct Line {
	Point start;
	Point end;

	inline Line() {}
	inline Line(const Point &s, const Point &e) :
			start(s),
			end(e) {}
} Line;

typedef struct Ray {
	Point origin;
	vec3 direction;

	inline Ray() :
			direction(0.0f, 0.0f, 1.0f) {}
	inline Ray(const Point &o, const vec3 &d) :
			origin(o),
			direction(d) {
		NormalizeDirection();
	}
	inline void NormalizeDirection() {
		Normalize(direction);
	}
} Ray;

typedef struct Sphere {
	Point position;
	float radius;

	inline Sphere() :
			radius(1.0f) {}
	inline Sphere(const Point &p, float r) :
			position(p),
			radius(r) {}
} Sphere;

typedef struct AABB {
	Point position;
	vec3 size;

	inline AABB() :
			size(1, 1, 1) {}
	inline AABB(const Point &p, const vec3 &s) :
			position(p),
			size(s) {}
} AABB;

typedef struct OBB {
	Point position;
	vec3 size;

	mat3 orientation;

	inline OBB() :
			size(1, 1, 1) {}
	inline OBB(const Point &p, const vec3 &s) :
			position(p),
			size(s) {}
	inline OBB(const Point &p, const vec3 &s, const mat3 &o) :
			position(p),
			size(s),
			orientation(o) {}
} OBB;

typedef struct Plane {
	vec3 normal;
	float distance;

	inline Plane() :
			normal(1, 0, 0) {}
	inline Plane(const vec3 &n, float d) :
			normal(n),
			distance(d) {}
} Plane;

typedef struct Triangle {
	union {
		struct {
			Point a;
			Point b;
			Point c;
		};
	};
	Point points[3];
	float values[9];

	inline Triangle() {}
	inline Triangle(const Point &p1, const Point &p2,
			const Point &p3) :
			a(p1),
			b(p2),
			c(p3) {}
};

typedef struct Interval {
	float min;
	float max;
} Interval;

typedef struct Mesh {
	int numTriangles;
	union {
		Triangle *triangles; //size = numTriangles
		Point *vertices; // size = numTriangles * 3
		float *values; //size = numTriangles *3 * 3
	};
	BVHNode *accelerator;
	Mesh() :
			numTriangles(0),
			values(0),
			accelerator(0) {}
} Mesh;

typedef struct BVHNode {
	AABB bounds;
	BVHNode *children;
	int numTriangles;
	int *triangles;
	BVHNode() :
			children(0),
			numTriangles(0),
			triangles(0) {}
} BVHNode;

typedef struct Frustum {
	union {
		struct {
			Plane top;
			Plane bottom;
			Plane left;
			Plane right;
			Plane near;
			Plane far;
		};
		Plane planes[6];
	};
	inline Frustum() {}
} Frustum;

typedef struct RaycastResult {
	vec3 point;
	vec3 normal;
	float t;
	bool hit;
} RaycastResult;
void ResetRaycastResult(RaycastResult *outResult);

typedef struct CollisionManifold {
	bool colliding;
	vec3 normal;
	float depth;
	std::vector<vec3> contacts;
} CollisionManifold;

void ResetCollisionManifold(CollisionManifold *result);

CollisionManifold FindCollisionFeatures(const Sphere &A,
		const Sphere &B);

CollisionManifold FindCollisionFeatures(const OBB &A,
		const Sphere &B);

std::vector<Point> GetVertices(const OBB &obb);
std::vector<Line> GetPlanes(const OBB &obb);
std::vector<Plane> GetPlanes(const OBB& obb);
bool ClipToPlane(const Plane& plane,
	const Line& line, Point* outPoint);
std::vector<Point> ClipEdgesToOBB(
		const std::vector<Line> &edges, const OBB &obb);
float PenetrationDepth(const OBB &o1, const OBB &o2,
		const vec3 &axis, bool *outShouldFlip);

CollisionManifold FindCollisionFeatures(const OBB &A,
		const OBB &B);

		float Length(const Line &line);
float LengthSq(const Line &line);

Ray FromPoints(const Point &from, const Point &to);

vec3 GetMin(const AABB &aabb);
vec3 GetMax(const AABB &aabb);
AABB FromMinMax(const vec3 &min, const vec3 &max);

float PlaneEquation(const Point &pt, const Plane &plane);

bool PointInSphere(const Point &point,
		const Sphere &sphere);
Point ClosestPoint(const Sphere &sphere,
		const Point &point);

bool PointInAABB(const Point &point, const AABB &aabb);
Point ClosestPoint(const AABB &aabb, const Point &point);

bool PointInOBB(const Point &point, const OBB &obb);
Point ClosestPoint(const OBB &obb, const Point &point);

bool PointOnPlane(const Point &point, const Plane &plane);
Point ClosestPoint(const Plane &plane, const Point &point);

bool PointOnLine(const Point &point, const Line &line);
Point ClosestPoint(const Line &line, const Point &point);

bool PointOnRay(const Point &point, const Ray &ray);
Point ClosestPoint(const Ray &ray, const Point &point);

bool SphereSphere(const Sphere &s1, const Sphere &s2);

bool SphereAABB(const Sphere &sphere, const AABB &aabb);

bool SphereOBB(const Sphere &sphere, const OBB &obb);

bool SpherePlane(const Sphere &sphere, const Plane &plane);

bool AABBAABB(const AABB &aabb1, const AABB &aabb2);

bool AABBOBB(const AABB &aabb, const OBB &obb);

Interval GetInterval(const AABB &rect, const vec3 &axis);
Interval GetInterval(const OBB &rect, const vec3 &axis);
bool OverlapOnAxis(const AABB &aabb, const OBB &obb,
		const vec3 &axis);

bool AABBPlane(const AABB &aabb, const Plane &plane);

bool OverlapOnAxis(const OBB &obb1, const OBB &obb2,
		const vec3 &axis);
bool OBBOBB(const OBB &obb1, const OBB &obb2);

bool OBBPlane(const OBB &obb, const Plane &plane);

bool PlanePlane(const Plane &plane1, const Plane &plane2);

bool Raycast(const Sphere &sphere,
		const Ray &ray, RaycastResult *outResult);

bool Raycast(const AABB &aabb, const Ray &ray,
		RaycastResult *outResult);

bool Raycast(const OBB &obb, const Ray &ray,
		RaycastResult *outResult);

bool Raycast(const Plane &plane, const Ray &ray,
		RaycastResult *outResult);

bool Linetest(const Sphere &sphere, const Line &line);

bool Linetest(const AABB &aabb, const Line &line);

bool Linetest(const OBB &obb, const Line &line);

bool Linetest(const Plane &plane, const Line &line);

bool PointInTriangle(const Point &p, const Triangle &t);

Plane FromTriangle(const Triangle &t);
Point ClosestPoint(const Triangle &t, const Point &p);

bool TriangleSphere(const Triangle &t, const Sphere &s);

Interval GetInterval(const Triangle &triangle, vec3 &axis);
bool OverlapOnAxis(const AABB &aabb, const Triangle &triangle,
		vec3 &axis);
bool TriangleAABB(const Triangle &t, const AABB &a);

bool OverlapOnAxis(const OBB &obb,
		const Triangle &triangle, const vec3 &axis);
bool TriangleOBB(const Triangle &t, const OBB &o);

bool TrianglePlane(const Triangle &t, const Plane &p);

bool OverlapOnAxis(const Triangle &t1,
		const Triangle &t2, const vec3 &axis);

bool TriangleTriangle(const Triangle &t1,
		const Triangle &t2);

//A-Edge /Triangle 0, point 0
//B-Edge / Triangle 0, point 1
//C-Edge / Triangle 1, point 0
//D-Edge / Triangle 1, point 1
vec3 SatCrossEdge(const vec3 &a, const vec3 &b,
		const vec3 &c, const vec3 &d);

bool TriangleTriangleRobust(const Triangle &t1,
		const Triangle &t2);

vec3 Barycentric(const Point &p, const Triangle &t);

float Raycast(const Triangle &triangle, const Ray &ray,
		RaycastResult *outResult);

bool Linetest(const Triangle &triangle, const Line &line);

void AccelerateMesh(Mesh &mesh);

void SplitBVHNode(BVHNode *node, const Mesh &mode,
		int depth);

void FreeBVHNode(BVHNode *node);

float MeshRay(const Mesh &mesh, const Ray &ray);
bool MeshAAB(const Mesh &mesh, const AABB &aabb);
//additional tests included with doanladable source code to investigate
float MeshRay(const Mesh &mesh, const Ray &ray);

bool Linetest(const Mesh &mesh, const Line &line);

bool MeshSphere(const Mesh &mesh, const Sphere &sphere);

bool MeshOBB(const Mesh &mesh, const OBB &obb);

bool MeshPlane(const Mesh &mesh, const Plane &plane);

bool MeshTriangle(const Mesh &mesh, const Triangle &triangle);

mat4 GetWorldMatrix(const Model &model);
OBB GetOBB(const Model &model);

float ModelRay(const Model &model, const Ray &ray);
bool Linetest(const Model &model, const Line &line);
bool ModelSphere(const Model &model, const Sphere &sphere);
bool ModelAABB(const Model &model, const AABB &aabb);
bool ModelOBB(const Model &model, const OBB &obb);
bool ModelPlane(const Model &model, const Plane &plane);
bool ModelTriangle(const Model &model,
		const Triangle &triangle);

Point Intersection(Plane p1, Plane p2, Plane p3);
void GetCorners(const Frustum &f, vec3 *outCorners);

bool Intersects(const Frustum &f, const Point &p);
bool Intersects(const Frustum &f, const Sphere &s);

float Classify(const AABB &aabb, const Plane &plane);
float Classify(const OBB &obb, const Plane &plane);
bool Intersects(const Frustum &f, const AABB &aabb);
bool Intersects(const Frustum &f, const OBB &obb);

vec3 Unproject(const vec3 &viewPoint,
		const vec2 &viewportOrigin, const vec2 &viewportSize,
		const mat4 &view, const mat4 &projection);
Ray GetPickRay(const vec2 &viewportPoint,
		const vec2 &viewportOrigin, const vec2 &viewportSize,
		const mat4 &view, const mat4 &projection);

class Model {
protected:
	Mesh *content;
	AABB bounds;

public:
	vec3 position;
	vec3 rotation;
	Model *parent;

	inline Model() :
			parent(0),
			content(0) {}
	inline Mesh *GetMesh() const {
		return content;
	}

	inline AABB GetBounds() const {
		return bounds;
	}
	void SetContent(Mesh *mesh);
};
#endif
