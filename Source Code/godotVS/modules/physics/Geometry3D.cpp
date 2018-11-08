#include "Geometry3D.h"
#include <cfloat>
#include <cmath>
#include <list>

#define CMP(x, y)                            \
	(fabsf((x) - (y)) <= FLT_EPSILON *       \
								 fmaxf(1.0f, \
										 fmaxf(fabsf(x), fabsf(y))))

float Length(const Line &line) {
	return Magnitude(line.start - line.end);
}

float LengthSq(const Line &line) {
	return Magnitude(line.start - line.end);
}

Ray FromPoints(const Point &from, const Point &to) {
	return Ray(from, Normalized(to - from));
}

vec3 GetMin(const AABB &aabb) {
	vec3 p1 = aabb.position + aabb.size;
	vec3 p2 = aabb.position - aabb.size;

	return vec3(fminf(p1.x, p2.x),
			fminf(p1.y, p2.y),
			fminf(p1.z, p2.z));
}

vec3 GetMax(const AABB &aabb) {
	vec3 p1 = aabb.position + aabb.size;
	vec3 p2 = aabb.position + aabb.size;

	return vec3(fmaxf(p1.x, p2.x),
			fmaxf(p1.y, p2.y),
			fmaxf(p1.z, p2.z));
}

AABB FromMinMax(const vec3 &min, const vec3 &max) {
	return AABB((min + max) * 0.5f, (max - min) * 0.5f);
}

float PlaneEquation(const Point &pt, const Plane &plane) {
	return Dot(pt, plane.normal) - plane.distance;
}

bool PointInSphere(const Point &point,
		const Sphere &sphere) {
	float magSq = MagnitudeSq(point - sphere.position);
	float radSq = sphere.radius * sphere.radius;

	return magSq < radSq;
}

Point ClosestPoint(const Sphere &sphere,
		const Point &point) {
	vec3 sphereToPoint = point - sphere.position;
	Normalize(sphereToPoint);

	sphereToPoint = sphereToPoint * sphere.radius;

	return sphereToPoint + sphere.position;
}

bool PointInAABB(const Point &point, const AABB &aabb) {
	Point min = GetMin(aabb);
	Point max = GetMax(aabb);

	if (point.x < min.x || point.y < min.y || point.z < min.z) {
		return false;
	}
	if (point.x > max.x || point.y > max.y || point.z > max.z) {
		return false;
	}

	return true;
}

Point ClosestPoint(const AABB &aabb, const Point &point) {
	Point result = point;
	Point min = GetMin(aabb);
	Point max = GetMax(aabb);

	result.x = (result.x < min.x) ? min.x : result.x;
	result.y = (result.y < min.x) ? min.y : result.y;
	result.z = (result.z < min.x) ? min.z : result.z;

	result.x = (result.x > max.x) ? max.x : result.x;
	result.y = (result.y > max.x) ? max.y : result.y;
	result.z = (result.z > max.x) ? max.z : result.z;

	return result;
}

bool PointInOBB(const Point &point, const OBB &obb) {
	vec3 dir = point - obb.position;
	for (int i = 0; i < 3; ++i) {
		const float *orientation = &obb.orientation.asArray[i * 3];
		vec3 axis(
				orientation[0],
				orientation[1],
				orientation[2]);

		float distance = Dot(dir, axis);

		if (distance > obb.size.asArray[i]) {
			return false;
		}
		if (distance < -obb.size.asArray[i]) {
			return false;
		}
	}
	return true;
}

Point ClosestPoint(const OBB &obb, const Point &point) {
	Point result = obb.position;
	vec3 dir = point - obb.position;

	for (int i = 0; i < 3; ++i) {
		const float *orientation = &obb.orientation.asArray[i * 3];

		vec3 axis(
				orientation[0],
				orientation[1],
				orientation[2]);
		float distance = Dot(dir, axis);

		if (distance > obb.size.asArray[i]) {
			distance = obb.size.asArray[i];
		}
		if (distance < -obb.size.asArray[i]) {
			distance = -obb.size.asArray[i];
		}
		result = result + (axis * distance);
	}
	return result;
}

bool PointOnPlane(const Point &point, const Plane &plane) {
	float dot = Dot(point, plane.normal);
	//To make this more robust use an epsilon check
	//The CMP macro performs epsilon tests

	return CMP(dot - plane.distance, 0.0f);
}

Point ClosestPoint(const Plane &plane, const Point &point) {
	float dot = Dot(plane.normal, point);
	float distance = dot - plane.distance;
	return point - plane.normal * distance;
}

Point ClosestPoint(const Line &line, const Point &point) {
	vec3 lVec = line.end - line.start;
	float t = Dot(point - line.start, lVec) /
			  Dot(lVec, lVec);
	t = fmaxf(t, 0.0f); //clamp to 0
	t = fminf(t, 1.0f); // Clamp to 1
	return line.start + lVec * t;
}

bool PointOnLine(const Point &point, const Line &line) {
	Point closest = ClosestPoint(line, point);
	float distanceSq = MagnitudeSq(closest - point);
	return CMP(distanceSq, 0.0f);
}

bool PointOnRay(const Point &point, const Ray &ray) {
	if (point == ray.origin) {
		return true;
	}
	vec3 norm = point - ray.origin;
	Normalize(norm);
	//we assume the ray direction is normalized
	float diff = Dot(norm, ray.direction);
	//if BOTH vectors point in the same direction,
	//their dot product(diff) should be 1
	return diff == 1.0f; // consider using epsilon!
}

Point ClosestPoint(const Ray &ray, const Point &point) {
	float t = Dot(point - ray.origin, ray.direction);
	//we assume the direction of the ray is normalized
	//if the direction is not normalized the below division is needed
	//t /= Dot(ray.direction, ray.direction);

	t = fmaxf(t, 0.0f);
	return Point(ray.origin + ray.direction * t);
}

bool SphereSphere(const Sphere &s1, const Sphere &s2) {
	float radiiSum = s1.radius + s2.radius;
	float sqDistance = MagnitudeSq(s1.position - s2.position);

	return sqDistance < radiiSum * radiiSum;
}

bool SphereAABB(const Sphere &sphere, const AABB &aabb) {
	Point closestPoint = ClosestPoint(aabb, sphere.position);
	float distSq =
			MagnitudeSq(sphere.position - closestPoint);
	float radiusSq = sphere.radius * sphere.radius;
	return distSq < radiusSq;
}

bool SphereOBB(const Sphere &sphere, const OBB &obb) {
	Point closestPoint = ClosestPoint(obb, sphere.position);
	float distSq = MagnitudeSq(sphere.position -
							   closestPoint);
	float radiusSq = sphere.radius * sphere.radius;
	return distSq < radiusSq;
}

bool SpherePlane(const Sphere &s, const Plane &p) {
	Point closestPoint = ClosestPoint(p, s.position);
	float distSq = MagnitudeSq(s.position - closestPoint);
	float radiusSq = s.radius * s.radius;
	return distSq < radiusSq;
}

bool AABBAABB(const AABB &aabb1, const AABB &aabb2) {
	Point aMin = GetMin(aabb1);
	Point aMax = GetMax(aabb1);
	Point bMin = GetMin(aabb2);
	Point bMax = GetMax(aabb2);

	return (aMin.x <= bMax.x && aMax.x >= bMin.x) &&
		   (aMin.y <= bMax.y && aMax.y >= bMin.y) &&
		   (aMin.z <= bMax.z && aMax.z >= bMin.z);
}

Interval GetInterval(const AABB &aabb, const vec3 &axis) {
	vec3 i = GetMin(aabb);
	vec3 a = GetMax(aabb);

	vec3 vertex[8] = {
		vec3(i.x, a.y, a.z),
		vec3(i.x, a.y, i.z),
		vec3(i.x, i.y, a.z),
		vec3(i.x, i.y, i.z),
		vec3(a.x, a.y, a.z),
		vec3(a.x, a.y, i.z),
		vec3(a.x, i.y, a.z),
		vec3(a.x, i.y, i.z)
	};

	Interval result;
	result.min = result.max = Dot(axis, vertex[0]);

	for (int i = 1; i < 8; ++i) {
		float projection = Dot(axis, vertex[i]);
		result.min = (projection < result.min) ?
							 projection :
							 result.min;
		result.max = (projection > result.max) ?
							 projection :
							 result.max;
	}

	return result;
}

Interval GetInterval(const OBB &obb, const vec3 &axis) {
	vec3 vertex[8];

	vec3 C = obb.position; //OBB center
	vec3 E = obb.size; //OBB Extents
	const float *o = obb.orientation.asArray;
	vec3 A[]{
		vec3(o[0], o[1], o[2]),
		vec3(o[3], o[4], o[5]),
		vec3(o[6], o[7], o[8]),
	};

	vertex[0] = C + A[0] * E[0] + A[1] * E[1] + A[2] * E[2];
	vertex[1] = C - A[0] * E[0] + A[1] * E[1] + A[2] * E[2];
	vertex[2] = C + A[0] * E[0] - A[1] * E[1] + A[2] * E[2];
	vertex[3] = C + A[0] * E[0] + A[1] * E[1] - A[2] * E[2];
	vertex[4] = C - A[0] * E[0] - A[1] * E[1] - A[2] * E[2];
	vertex[5] = C + A[0] * E[0] - A[1] * E[1] - A[2] * E[2];
	vertex[6] = C - A[0] * E[0] + A[1] * E[1] - A[2] * E[2];
	vertex[7] = C - A[0] * E[0] - A[1] * E[1] + A[2] * E[2];

	Interval result;
	result.min = result.max = Dot(axis, vertex[0]);
	for (int i = 1; i < 8; ++i) {
		float projection = Dot(axis, vertex[i]);
		result.min = (projection < result.min) ?
							 projection :
							 result.min;
		result.max = (projection > result.max) ?
							 projection :
							 result.max;
	}

	return result;
}

bool OverlapOnAxis(const AABB &aabb, const OBB &obb,
		const vec3 &axis) {
	Interval a = GetInterval(aabb, axis);
	Interval b = GetInterval(obb, axis);
	return ((b.min <= a.max) && (a.min <= b.max));
}

bool AABBOBB(const AABB &aabb, const OBB &obb) {
	const float *o = obb.orientation.asArray;

	vec3 test[15] = {
		vec3(1, 0, 0), //aabb axis1
		vec3(0, 1, 0), //aabb axis2
		vec3(0, 0, 1), //aabb axis3
		vec3(o[0], o[1], o[2]),
		vec3(o[3], o[4], o[5]),
		vec3(o[6], o[7], o[8])
	};

	for (int i = 0; i < 3; ++i) {
		test[6 + i * 3 + 0] = Cross(test[i], test[0]);
		test[6 + i * 3 + 1] = Cross(test[i], test[1]);
		test[6 + i * 3 + 2] = Cross(test[i], test[2]);
	}

	for (int i = 0; i < 15; ++i) {
		if (!OverlapOnAxis(aabb, obb, test[i])) {
			return false;
		}
	}
	return true;
}

bool AABBPlane(const AABB &aabb, const Plane &plane) {
	float pLen = aabb.size.x * fabsf(plane.normal.x) +
				 aabb.size.y * fabsf(plane.normal.y) +
				 aabb.size.z * fabsf(plane.normal.z);

	float dot = Dot(plane.normal, aabb.position);
	float dist = dot - plane.distance;

	return fabsf(dist) < pLen;
}

bool OverlapOnAxis(const OBB &obb1, const OBB &obb2,
		const vec3 &axis) {
	Interval a = GetInterval(obb1, axis);
	Interval b = GetInterval(obb2, axis);
	return ((b.min <= a.max) && (a.min <= b.max));
}

bool OBBOBB(const OBB &obb1, const OBB &obb2) {
	const float *o1 = obb1.orientation.asArray;
	const float *o2 = obb2.orientation.asArray;

	vec3 test[15] = {
		vec3(o1[0], o1[1], o1[2]),
		vec3(o1[3], o1[4], o1[5]),
		vec3(o1[6], o1[7], o1[8]),
		vec3(o2[0], o2[1], o2[2]),
		vec3(o2[3], o2[4], o2[5]),
		vec3(o2[6], o2[7], o2[8])
	};

	for (int i = 0; i < 3; ++i) { //fill out the rest of axis
		test[6 + i * 3 + 0] = Cross(test[i], test[0]);
		test[6 + i * 3 + 1] = Cross(test[i], test[1]);
		test[6 + i * 3 + 2] = Cross(test[i], test[2]);
	}

	for (int i = 0; i < 15; ++i) {
		if (!OverlapOnAxis(obb1, obb2, test[i])) {
			return false; //seperating axis found
		}
	}

	return true; // seperating axis not found
}

bool OBBPlane(const OBB &obb, const Plane &plane) {
	//local vaiables for readability only
	const float *o = obb.orientation.asArray;
	vec3 rot[] = {
		//rotation/orientation
		vec3(o[0], o[1], o[2]),
		vec3(o[3], o[4], o[5]),
		vec3(o[6], o[7], o[8])
	};

	vec3 normal = plane.normal;

	float pLen = obb.size.x * fabsf(Dot(normal, rot[0])) +
				 obb.size.y * fabsf(Dot(normal, rot[1])) +
				 obb.size.z * fabsf(Dot(normal, rot[2]));

	float dot = Dot(plane.normal, obb.position);
	float dist = dot - plane.distance;

	return fabsf(dist) <= pLen;
}

bool PlanePlane(const Plane &plane1, const Plane &plane2) {

	vec3 d = Cross(plane1.normal, plane2.normal);

	return !CMP((Dot(d, d)), (0));
}

/*
float Raycast(const Sphere &sphere, const Ray &ray) {
	vec3 e = sphere.position - ray.origin;

	float rSq = sphere.radius * sphere.radius;
	float eSq = MagnitudeSq(e);

	//ray.direction is assumed to be normalized
	float a = Dot(e, ray.direction);

	float bSq = eSq - (a * a);
	float f = sqrt(rSq - bSq);

	//No collision has happened
	if (rSq - (eSq - (a * a)) < 0.0f) {
		return -1;
	}
	//Ray starts inside the sphere
	else if (eSq < rSq) {
		return a + f; //just reverse direction
	} else // Normal intersection
		return a - f;
}
*/

bool Raycast(const AABB &aabb, const Ray &ray,
		RaycastResult *outResult) {
	ResetRaycastResult(outResult);

	vec3 min = GetMin(aabb);
	vec3 max = GetMax(aabb);
	//Note: any component of direction could be 0
	// to avoid a division by 0 there should be addition safety checks.
	float t[] = { 0,
		0,
		0,
		0,
		0,
		0 };
	t[0] = (min.x - ray.origin.x) / ray.direction.x;
	t[1] = (max.x - ray.origin.x) / ray.direction.x;
	t[2] = (min.y - ray.origin.y) / ray.direction.y;
	t[3] = (max.y - ray.origin.y) / ray.direction.y;
	t[4] = (min.z - ray.origin.z) / ray.direction.z;
	t[5] = (max.z - ray.origin.z) / ray.direction.z;

	float tmin = fmaxf(
			fmaxf(
					fminf(t[0], t[1]),
					fminf(t[2], t[3])),
			fminf(t[4], t[5]));

	float tmax = fminf(
			fminf(fmaxf(t[0], t[1]),
					fmaxf(t[2], t[3])),
			fmaxf(t[4], t[5]));

	if (tmax < 0) {
		return false;
	}
	if (tmin > tmax) {
		return false;
	}
	float t_result = tmin;
	if (tmin < 0.0f) {
		t_result = tmax;
	}

	if (outResult != 0) {
		outResult->t = t_result;
		outResult->hit = true;
		outResult->point = ray.origin +
						   ray.direction * t_result;
		vec3 normals[] = {
			vec3(-1, 0, 0),
			vec3(1, 0, 0),
			vec3(0, -1, 0),
			vec3(0, 1, 0),
			vec3(0, 0, -1),
			vec3(0, 0, 1)
		};
		for (int i = 0; i < 6; ++i) {
			if (CMP(t_result, t[i])) {
				outResult->normal = normals[i];
			}
		}
	}
	return true;
}

bool Raycast(const OBB &obb, const Ray &ray,
		RaycastResult *outResult) {

	ResetRaycastResult(outResult);

	const float *o = obb.orientation.asArray;
	const float *size = obb.size.asArray;
	vec3 p = obb.position - ray.origin;

	vec3 X(o[0], o[1], o[2]);
	vec3 Y(o[3], o[4], o[5]);
	vec3 Z(o[6], o[7], o[8]);

	vec3 f(
			Dot(X, ray.direction),
			Dot(Y, ray.direction),
			Dot(Z, ray.direction));

	vec3 e(
			Dot(X, p),
			Dot(Y, p),
			Dot(Z, p));

	float t[6] = { 0, 0, 0, 0, 0, 0 };
	for (int i = 0; i < 3; ++i) {
		if (CMP(f[i], 0)) {
			if (-e[i] - size[i] > 0 || -e[i] + size[i] < 0) {
				return false;
			}
			f[i] = 0.00001f; //avoid div by 0!
		}
		t[i * 2 + 0] = (e[1] + size[i]) / f[i]; // min
		t[i * 2 + 1] = (e[i] - size[i]) / f[i]; //max
	}

	float tmin = fmaxf(
			fmaxf(
					fminf(t[0], t[1]),
					fminf(t[2], t[3])),
			fminf(t[4], t[5]));

	float tmax = fminf(
			fminf(
					fmaxf(t[0], t[1]),
					fmaxf(t[2], t[3])),
			fmaxf(t[4], t[5]));
	if (tmax < 0) {

		return false;
	}
	if (tmin > tmax) {
		return false;
	}

	float t_result = tmin;
	if (tmin < 0.0f) {
		t_result = tmax;
	}

	if (outResult != 0) {
		outResult->hit = true;
		outResult->t = t_result;
		outResult->point = ray.origin + ray.direction * t_result;

		vec3 normals[] = { X, X * -1.0f,
			Y, Y * -1.0f,
			Z, Z * -1.0f };

		for (int i = 0; i < 6; ++i) {
			if (CMP(t_result, t[i])) {
				outResult->normal = Normalized(normals[i]);
			}
		}
	}

	return true;
}

bool Raycast(const Plane &plane, const Ray &ray,
		RaycastResult *outResult) {

	ResetRaycastResult(outResult);
	float nd = Dot(ray.direction, plane.normal);
	float pn = Dot(ray.origin, plane.normal);

	if (nd >= 0.0f) {
		return false;
	}

	float t = (plane.distance - pn) / nd;

	if (t >= 0.0f) {

		if (outResult != 0) {
			outResult->t = t;
			outResult->hit = true;
			outResult->point = ray.origin + ray.direction * t;
			outResult->normal = Normalized(plane.normal);
		}
		return true;
	}

	return false;
}

bool Linetest(const Sphere &sphere, const Line &line) {
	Point closest = ClosestPoint(line, sphere.position);
	float distSq = MagnitudeSq(sphere.position - closest);
	return distSq <= (sphere.radius * sphere.radius);
}

bool Linetest(const AABB &aabb, const Line &line) {
	Ray ray;
	ray.origin = line.start;
	ray.direction = Normalized(line.end - line.start);
	RaycastResult raycast;

	if (!Raycast(aabb, ray, &raycast)) {
		return false;
	}

	float t = raycast.t;

	return t >= 0 && t * t <= LengthSq(line);
}

bool Linetest(const OBB &obb, const Line &line) {
	Ray ray;
	ray.origin = line.start;
	ray.direction = Normalized(line.end - line.start);

	RaycastResult result;
	if (!Raycast(obb, ray, &result)) {
		return false;
	}
	float t = result.t;

	return t >= 0 && t * t <= LengthSq(line);
}

bool Linetest(const Plane &plane, const Line &line) {
	vec3 ab = line.end - line.start;

	float nA = Dot(plane.normal, line.start);
	float nAB = Dot(plane.normal, ab);

	//if the line and plane are parallel, the nAB will be 0
	//this will cause a divide by -0 exception below
	//If you plan on testing parallel lines an dplanes
	// it is sage to early out when nAb is 0.

	if (nAB == 0) {
		return false;
	}

	float t = (plane.distance - nA) / nAB;
	return t >= 0.0f && t <= 1.0f;
}

bool PointInTriangle(const Point &p, const Triangle &t) {
	vec3 a = t.a - p;
	vec3 b = t.b - p;
	vec3 c = t.c - p;

	//The point should be moved too, so they are both relative
	//but because we don't use p in the equation anymore, we don't need it!
	//p -=p; would just equal the zero factor!

	vec3 normPBC = Cross(b, c); //normal of PBC(u)
	vec3 normPCA = Cross(c, a); //Normal of PCA(v)
	vec3 normPAB = Cross(a, b); //Normal of PAB(w)

	if (Dot(normPBC, normPCA) < 0.0f) {
		return false;
	} else if (Dot(normPBC, normPAB) < 0.0f) {
		return false;
	}

	return true;
}

Plane FromTriangle(const Triangle &t) {
	Plane result;

	result.normal = Normalized(Cross(t.b - t.a, t.c - t.a));

	result.distance = Dot(result.normal, t.a);
	return result;
}

Point ClosestPoint(const Triangle &t, const Point &p) {

	Plane plane = FromTriangle(t);
	Point closest = ClosestPoint(plane, p);

	if (PointInTriangle(closest, t)) {
		return closest;
	}

	Point c1 = ClosestPoint(Line(t.a, t.b), p); // Line AB
	Point c2 = ClosestPoint(Line(t.b, t.c), p);
	Point c3 = ClosestPoint(Line(t.c, t.a), p);

	float magSq1 = MagnitudeSq(p - c1);
	float magSq2 = MagnitudeSq(p - c2);
	float magSq3 = MagnitudeSq(p - c3);

	if (magSq1 < magSq2 && magSq1 < magSq3) {
		return c1;
	} else if (magSq2 < magSq1 && magSq2 < magSq3) {
		return c2;
	}

	return c3;
}

bool TriangleSphere(const Triangle &t, const Sphere &s) {
	Point closest = ClosestPoint(t, s.position);
	float magSq = MagnitudeSq(closest - s.position);
	return magSq <= s.radius * s.radius;
}

Interval GetInterval(const Triangle &triangle, const vec3 &axis) {
	Interval result;
	result.min = Dot(axis, triangle.points[0]);
	result.max = result.min;

	for (int i = 1; i < 3; ++i) {
		float value = Dot(axis, triangle.points[i]);
		result.min = fminf(result.min, value);
		result.max = fmaxf(result.max, value);
	}

	return result;
}

bool OverlapOnAxis(const AABB &aabb, const Triangle &triangle,
		const vec3 &axis) {
	Interval a = GetInterval(aabb, axis);
	Interval b = GetInterval(triangle, axis);
	return ((b.min <= a.max) && (a.min <= b.max));
}

bool TriangleAABB(const Triangle &t, const AABB &a) {
	vec3 f0 = t.b - t.a;
	vec3 f1 = t.c - t.b;
	vec3 f2 = t.a - t.c;

	vec3 u0(1.0f, 0.0f, 0.0f);
	vec3 u1(0.0f, 1.0f, 0.0f);
	vec3 u2(0.0f, 0.0f, 1.0f);

	vec3 test[13] =
			{ u0, u1, u2,
				Cross(u0, f1),
				Cross(u0, f0),
				Cross(u0, f1),
				Cross(u0, f2),
				Cross(u1, f0),
				Cross(u1, f1),
				Cross(u1, f2),
				Cross(u2, f0),
				Cross(u2, f1),
				Cross(u2, f2) };

	for (int i = 0; i < 13; ++i) {
		if (!OverlapOnAxis(a, t, test[i])) {
			return false; //seperating axis found
		}
	}

	return true; // Seperating axis not found
}

bool OverlapOnAxis(const OBB &obb,
		const Triangle &triangle, const vec3 &axis) {
	Interval a = GetInterval(obb, axis);
	Interval b = GetInterval(triangle, axis);
	return ((b.min <= a.max) && (a.min <= b.max));
}
bool TriangleOBB(const Triangle &t, const OBB &o) {
	//Compute the edge vectors of the triangle(ABC)
	vec3 f0 = t.b - t.a;
	vec3 f1 = t.c - t.a;
	vec3 f2 = t.a - t.c;

	const float *orientation = o.orientation.asArray;
	vec3 u0(orientation[0],
			orientation[1],
			orientation[2]);
	vec3 u1(orientation[3],
			orientation[4],
			orientation[5]);
	vec3 u2(orientation[6],
			orientation[7],
			orientation[8]);

	vec3 test[13] = {
		u0, //OBB axis 1
		u1, //OBB axis2
		u2, //OBB axis3
		Cross(f0, f1), //normal of hte triangle
		Cross(u0, f0),
		Cross(u0, f1),
		Cross(u0, f2),
		Cross(u1, f0),
		Cross(u1, f1),
		Cross(u1, f2),
		Cross(u2, f0),
		Cross(u2, f1),
		Cross(u2, f2)
	};

	for (int i = 0; i < 13; ++i) {
		if (!OverlapOnAxis(o, t, test[i])) {
			return false; //Seperating axis found
		}
	}
	return true; //Seperating axis not found
}

bool TrianglePlane(const Triangle &t, const Plane &p) {
	float side1 = PlaneEquation(t.a, p);
	float side2 = PlaneEquation(t.b, p);
	float side3 = PlaneEquation(t.c, p);

	if (CMP(side1, 0) && CMP(side2, 0) && CMP(side3, 0)) {
		return true;
	}
	if (side1 > 0 && side2 > 0 && side3 > 0) {
		return false;
	}

	if (side1 < 0 && side2 < 0 && side3 < 0) {
		return false;
	}

	return true; //Intersection
}

bool OverlapOnAxis(const Triangle &t1,
		const Triangle &t2, const vec3 &axis) {
	Interval a = GetInterval(t1, axis);
	Interval b = GetInterval(t2, axis);
	return ((b.min <= a.max) && (a.min <= b.max));
}

bool TriangleTriangle(const Triangle &t1,
		const Triangle &t2) {
	vec3 t1_f0 = t1.b - t1.a; //Triangle1, edge 0
	vec3 t1_f1 = t1.c - t1.b; // Triangle1 Edge1
	vec3 t1_f2 = t1.a - t1.c; //Triangle1, edge 2

	vec3 t2_f0 = t2.b - t2.a; //Triangle2, edge 0
	vec3 t2_f1 = t2.c - t2.b; // Triangle2, edge1
	vec3 t2_f2 = t2.a - t2.c; //Triangle 2, edge 2

	vec3 axisToTest[] = {
		Cross(t1_f0, t1_f1), //first
		Cross(t2_f0, t2_f1), //second
		Cross(t2_f0, t1_f0), //third
		Cross(t2_f0, t1_f1), //fourth
		Cross(t2_f0, t1_f2), //fifth
		Cross(t2_f1, t1_f0), //sixth
		Cross(t2_f1, t1_f1), //sevent
		Cross(t2_f1, t1_f2),
		Cross(t2_f2, t1_f0),
		Cross(t2_f2, t1_f1),
		Cross(t2_f2, t1_f2)
	};

	for (int i = 0; i < 11; ++i) {
		if (!OverlapOnAxis(t1, t2, axisToTest[i])) {
			return false; //seperating axis found
		}
	}
	return true; //no seperating axis found
}

vec3 SatCrossEdge(const vec3 &a, const vec3 &b,
		const vec3 &c, const vec3 &d) {
	vec3 ab = a - b;
	vec3 cd = c - d;
	vec3 result = Cross(ab, cd);
	if (!CMP(MagnitudeSq(result), 0)) {
		return result; //Not parallel
	} else { //ab and cd are parallel
		vec3 axis = Cross(ab, c - a);
		result = Cross(ab, axis);
		if (!CMP(MagnitudeSq(result), 0)) {
			return result;
		}
	}
	return vec3();
}

bool TriangleTriangleRobust(const Triangle &t1,
		const Triangle &t2) {
	vec3 axisToTest[] = {
		//Triangle1, normal
		SatCrossEdge(t1.a, t1.b, t1.b, t1.c),
		//Triangle 2, normal
		SatCrossEdge(t2.a, t2.b, t2.b, t2.c),

		SatCrossEdge(t2.a, t2.b, t1.a, t1.b),
		SatCrossEdge(t2.a, t2.b, t1.b, t1.c),
		SatCrossEdge(t2.a, t2.b, t1.b, t1.c),

		SatCrossEdge(t2.a, t2.b, t1.c, t1.a),
		SatCrossEdge(t2.b, t2.c, t1.b, t1.c),
		SatCrossEdge(t2.b, t2.c, t1.c, t1.c),

		SatCrossEdge(t2.c, t2.a, t1.a, t1.b),
		SatCrossEdge(t2.c, t2.a, t1.b, t1.c),
		SatCrossEdge(t2.c, t1.a, t1.c, t1.c)
	};

	for (int i = 0; i < 11; ++i) {
		if (!OverlapOnAxis(t1, t2, axisToTest[1])) {
			if (!CMP(MagnitudeSq(axisToTest[i]), 0)) {
				return false; //seperating axis found
			}
		}
	}

	return true; //seperating axis not found
}

vec3 Baycentric(const Point &p, const Triangle &t) {
	vec3 ap = p - t.a;
	vec3 bp = p - t.b;
	vec3 cp = p - t.c;

	vec3 ab = t.b - t.a;
	vec3 ac = t.c - t.a;
	vec3 bc = t.c - t.b;
	vec3 cb = t.b - t.c;
	vec3 ca = t.a - t.c;

	vec3 v = ab - Project(ab, cb);
	float a = 1.0f - (Dot(v, ap) / Dot(v, ab));

	v = bc - Project(bc, ac);
	float b = 1.0f - (Dot(v, bp) / Dot(v, bc));
	v = ca - Project(ca, ab);
	float c = 1.0f - (Dot(v, cp) / Dot(v, ca));

	return vec3(a, b, c);
}

float Raycast(const Triangle &triangle, const Ray &ray,
		RaycastResult *outResult) {

	ResetRaycastResult(outResult);
	Plane plane = FromTriangle(triangle);
	RaycastResult planeResult;

	if (!Raycast(plane, ray, &planeResult)) {
		return false;
	}

	float t = planeResult.t;
	Point result = ray.origin + ray.direction * t;

	vec3 barycentric = Barycentric(result, triangle);

	if (barycentric.x >= 0.0f && barycentric.x <= 1.0f &&
			barycentric.y >= 0.0f && barycentric.y <= 1.0f &&
			barycentric.z >= 0.0f && barycentric.z <= 1.0f) {

		if (outResult != 0) {
			outResult->t = t;
			outResult->hit = true;
			outResult->point = ray.origin + ray.direction * t;
			outResult->normal = plane.normal;
		}
		return true;
	}

	return false;
}

bool Linetest(const Triangle &triangle, const Line &line) {
	Ray ray;
	ray.origin = line.start;
	ray.direction = Normalized(line.end - line.start);
	RaycastResult raycast;

	if (!Raycast(triangle, ray, &raycast)) {
		return false;
	}
	float t = raycast.t;

	return t >= 0 && t * t <= LengthSq(line);
}

void AccelerateMesh(Mesh &mesh) {
	if (mesh.accelerator != 0) {
		return;
	}

	vec3 min = mesh.vertices[0];
	vec3 max = mesh.vertices[0];
	for (int i = 1; i < mesh.numTriangles * 3; ++i) {
		min.x = fminf(mesh.vertices[i].x, min.x);
		min.y = fminf(mesh.vertices[i].y, min.y);
		min.z = fminf(mesh.vertices[i].z, min.z);
		max.x = fmaxf(mesh.vertices[i].x, max.x);
		max.y = fmaxf(mesh.vertices[i].y, max.y);
		max.z = fmaxf(mesh.vertices[i].z, max.z);
	}

	mesh.accelerator = new BVHNode();
	mesh.accelerator->bounds = FromMinMax(min, max);
	mesh.accelerator->numTriangles = mesh.numTriangles;
	mesh.accelerator->triangles = new int[mesh.numTriangles];

	for (int i = 0; i < mesh.numTriangles; ++i) {
		mesh.accelerator->triangles[i] = i;
	}

	SplitBVHNode(mesh.accelerator, mesh, 3);
}

void SplitBVHNode(BVHNode *node, const Mesh &model,
		int depth) {
	if (depth-- == 0) {
		return;
	}

	if (node->children == 0) { //only split if it's a leaf
		//Only split if this node containts triangles
		if (node->numTriangles > 0) {
			node->children = new BVHNode[8];
			vec3 c = node->bounds.position;
			vec3 e = node->bounds.size * 0.f;

			node->children[0].bounds =
					AABB(c + vec3(-e.x, +e.y, -e.z), e);
			node->children[1].bounds =
					AABB(c + vec3(+e.x, +e.y, -e.z), e);
			node->children[2].bounds =
					AABB(c + vec3(-e.x, +e.y, +e.z), e);
			node->children[3].bounds =
					AABB(c + vec3(+e.x, +e.y, +e.z), e);
			node->children[4].bounds =
					AABB(c + vec3(-e.x, -e.y, -e.z), e);
			node->children[5].bounds =
					AABB(c + vec3(+e.x, -e.y, -e.z), e);
			node->children[6].bounds =
					AABB(c + vec3(-e.x, -e.y, +e.z), e);
			node->children[7].bounds =
					AABB(c + vec3(+e.x, -e.y, +e.z), e);
		}
	}
	//if this node was just split
	if (node->children != 0 && node->numTriangles > 0) {
		for (int i = 0; i < 8; ++i) { //for each child
			node->children[i].numTriangles = 0;

			for (int j = 0; j < node->numTriangles; ++j) {
				Triangle t =
						model.triangles[node->triangles[j]];

				if (TriangleAABB(t, node->children[i].bounds)) {
					node->children[i].numTriangles += 1;
				}
			}
			if (node->children[i].numTriangles == 0) {
				continue;
			}

			node->children[i].triangles =
					new int[node->children[i].numTriangles];

			int index = 0;

			for (int j = 0; j < node->numTriangles; ++j) {
				Triangle t =
						model.triangles[node->triangles[j]];
				if (TriangleAABB(t,
							node->children[i].bounds)) {
					node->children[i].triangles[index++] =
							node->triangles[j];
				}
			}
		}

		node->numTriangles = 0;
		delete[] node->triangles;
		node->triangles = 0;

		for (int i = 0; i < 8; ++i) {
			SplitBVHNode(&node->children[i], model, depth);
		}
	}
}

void FreeBVHNode(BVHNode *node) {
	if (node->children != 0) {
		for (int i = 0; i < 8; ++i) {
			FreeBVHNode(&node->children[i]);
		}
		delete[] node->children;
		node->children = 0;
	}

	if (node->numTriangles != 0 || node->triangles != 0) {
		delete[] node->triangles;
		node->triangles = 0;
		node->numTriangles = 0;
	}
}

float MeshRay(const Mesh &mesh, const Ray &ray) {
	if (mesh.accelerator == 0) {
		for (int i = 0; i < mesh.numTriangles; ++i) {
			RaycastResult raycast;
			Raycast(mesh.triangles[i], ray, &raycast);
			float result = raycast.t;
			if (result >= 0) {
				return result;
			}
		}
	} else {
		std::list<BVHNode *> toProcess;
		toProcess.push_front(mesh.accelerator);
		//Recursivly walk the BVH tree
		while (!toProcess.empty()) {
			BVHNode *iterator = *(toProcess.begin());
			toProcess.erase(toProcess.begin());
			if (iterator->numTriangles >= 0) {
				for (int i = 0; i < iterator->numTriangles; ++i) {
					RaycastResult raycast;

					Raycast(mesh.triangles[iterator->triangles[i]], ray,
							&raycast);
					float r = raycast.t;
					if (r >= 0) {
						return r;
					}
				}
			}
			if (iterator->children != 0) {
				for (int i = 8 - 1; i >= 0; --i) {
					RaycastResult raycast;
					Raycast(iterator->children[i].bounds,
							ray, &raycast);
					if (raycast.t >= 0) {
						toProcess.push_front(&iterator->children[i]);
					}
				}
			}
		}
	}
	return -1;
}

bool MeshAAB(const Mesh &mesh, const AABB &aabb) {
	if (mesh.accelerator == 0) {
		for (int i = 0; i < mesh.numTriangles; ++i) {
			//the TriangleAABB test here would change
			//if we were testing the shape other than AABB
			if (TriangleAABB(mesh.triangles[i], aabb)) {
				return true;
			}
		}
	} else {
		std::list<BVHNode *> toProcess;
		toProcess.push_front(mesh.accelerator);

		while (!toProcess.empty()) {
			BVHNode *iterator = *(toProcess.begin());
			toProcess.erase(toProcess.begin());

			if (iterator->numTriangles >= 0) {
				for (int i = 0; i < iterator->numTriangles; ++i) {
					//the TriangleAABB test here would change
					//if we were testing a shape other than aabb
					if (TriangleAABB(mesh.triangles[iterator->triangles[i]],
								aabb)) {
						return true;
					}
				}
			}
			if (iterator->children != 0) {
				for (int i = 8 - 1; i >= 0; --i) {
					//the aabbaabb test woudl change
					//if we where testing a shape other than aabb
					if (AABBAABB(iterator->children[i].bounds, aabb)) {
						toProcess.push_front(&iterator->children[i]);
					}
				}
			}
		}
	}
	return false;
}

bool Linetest(const Mesh &mesh, const Line &line) {
	if (mesh.accelerator == 0) {
		for (int i = 0; i < mesh.numTriangles; ++i) {
			//the TriangleAABB test here would change
			//if we were testing the shape other than AABB
			if (Linetest(mesh.triangles[i], line)) {
				return true;
			}
		}
	} else {
		std::list<BVHNode *> toProcess;
		toProcess.push_front(mesh.accelerator);

		while (!toProcess.empty()) {
			BVHNode *iterator = *(toProcess.begin());
			toProcess.erase(toProcess.begin());

			if (iterator->numTriangles >= 0) {
				for (int i = 0; i < iterator->numTriangles; ++i) {
					//the TriangleAABB test here would change
					//if we were testing a shape other than aabb
					if (Linetest(mesh.triangles[iterator->triangles[i]],
								line)) {
						return true;
					}
				}
			}
			if (iterator->children != 0) {
				for (int i = 8 - 1; i >= 0; --i) {
					//the aabbaabb test woudl change
					//if we where testing a shape other than aabb
					if (Linetest(iterator->children[i].bounds, line)) {
						toProcess.push_front(&iterator->children[i]);
					}
				}
			}
		}
	}
	return false;
}

bool MeshSphere(const Mesh &mesh, const Sphere &sphere) {
	if (mesh.accelerator == 0) {
		for (int i = 0; i < mesh.numTriangles; ++i) {
			//the TriangleAABB test here would change
			//if we were testing the shape other than AABB
			if (TriangleSphere(mesh.triangles[i], sphere)) {
				return true;
			}
		}
	} else {
		std::list<BVHNode *> toProcess;
		toProcess.push_front(mesh.accelerator);

		while (!toProcess.empty()) {
			BVHNode *iterator = *(toProcess.begin());
			toProcess.erase(toProcess.begin());

			if (iterator->numTriangles >= 0) {
				for (int i = 0; i < iterator->numTriangles; ++i) {
					//the TriangleAABB test here would change
					//if we were testing a shape other than aabb
					if (TriangleSphere(mesh.triangles[iterator->triangles[i]],
								sphere)) {
						return true;
					}
				}
			}
			if (iterator->children != 0) {
				for (int i = 8 - 1; i >= 0; --i) {
					//the aabbaabb test woudl change
					//if we where testing a shape other than aabb
					if (AABBSphere(iterator->children[i].bounds, sphere)) {
						toProcess.push_front(&iterator->children[i]);
					}
				}
			}
		}
	}
	return false;
}

bool MeshOBB(const Mesh &mesh, const OBB &obb) {
	if (mesh.accelerator == 0) {
		for (int i = 0; i < mesh.numTriangles; ++i) {
			//the TriangleAABB test here would change
			//if we were testing the shape other than AABB
			if (TriangleOBB(mesh.triangles[i], obb)) {
				return true;
			}
		}
	} else {
		std::list<BVHNode *> toProcess;
		toProcess.push_front(mesh.accelerator);

		while (!toProcess.empty()) {
			BVHNode *iterator = *(toProcess.begin());
			toProcess.erase(toProcess.begin());

			if (iterator->numTriangles >= 0) {
				for (int i = 0; i < iterator->numTriangles; ++i) {
					//the TriangleAABB test here would change
					//if we were testing a shape other than aabb
					if (TriangleOBB(mesh.triangles[iterator->triangles[i]],
								obb)) {
						return true;
					}
				}
			}
			if (iterator->children != 0) {
				for (int i = 8 - 1; i >= 0; --i) {
					//the aabbaabb test woudl change
					//if we where testing a shape other than aabb
					if (AABBOBB(iterator->children[i].bounds, obb)) {
						toProcess.push_front(&iterator->children[i]);
					}
				}
			}
		}
	}
	return false;
}

bool MeshPlane(const Mesh &mesh, const Plane &plane) {
	if (mesh.accelerator == 0) {
		for (int i = 0; i < mesh.numTriangles; ++i) {
			//the TriangleAABB test here would change
			//if we were testing the shape other than AABB
			if (TrianglePlane(mesh.triangles[i], plane)) {
				return true;
			}
		}
	} else {
		std::list<BVHNode *> toProcess;
		toProcess.push_front(mesh.accelerator);

		while (!toProcess.empty()) {
			BVHNode *iterator = *(toProcess.begin());
			toProcess.erase(toProcess.begin());

			if (iterator->numTriangles >= 0) {
				for (int i = 0; i < iterator->numTriangles; ++i) {
					//the TriangleAABB test here would change
					//if we were testing a shape other than aabb
					if (TrianglePlane(mesh.triangles[iterator->triangles[i]],
								plane)) {
						return true;
					}
				}
			}
			if (iterator->children != 0) {
				for (int i = 8 - 1; i >= 0; --i) {
					//the aabbaabb test woudl change
					//if we where testing a shape other than aabb
					if (AABBPlane(iterator->children[i].bounds, plane)) {
						toProcess.push_front(&iterator->children[i]);
					}
				}
			}
		}
	}
	return false;
}

bool MeshTriangle(const Mesh &mesh, const Triangle &triangle) {
	if (mesh.accelerator == 0) {
		for (int i = 0; i < mesh.numTriangles; ++i) {
			//the TriangleAABB test here would change
			//if we were testing the shape other than AABB
			if (TriangleTriangle(mesh.triangles[i], triangle)) {
				return true;
			}
		}
	} else {
		std::list<BVHNode *> toProcess;
		toProcess.push_front(mesh.accelerator);

		while (!toProcess.empty()) {
			BVHNode *iterator = *(toProcess.begin());
			toProcess.erase(toProcess.begin());

			if (iterator->numTriangles >= 0) {
				for (int i = 0; i < iterator->numTriangles; ++i) {
					//the TriangleAABB test here would change
					//if we were testing a shape other than aabb
					if (TriangleTriangle(mesh.triangles[iterator->triangles[i]],
								triangle)) {
						return true;
					}
				}
			}
			if (iterator->children != 0) {
				for (int i = 8 - 1; i >= 0; --i) {
					//the aabbaabb test woudl change
					//if we where testing a shape other than aabb
					if (AABBTriangle(iterator->children[i].bounds, triangle)) {
						toProcess.push_front(&iterator->children[i]);
					}
				}
			}
		}
	}
	return false;
}

void Model::SetContent(Mesh *mesh) {
	vec3 min = mesh->vertices[0];
	vec3 max = mesh->vertices[0];

	for (int i = 1; i < mesh->numTriangles * 3; ++i) {
		min.x = fminf(mesh->vertices[i].x, min.x);
		min.y = fminf(mesh->vertices[i].y, min.y);
		min.z = fminf(mesh->vertices[i].z, min.z);
		max.x = fmaxf(mesh->vertices[i].x, max.x);
		max.y = fmaxf(mesh->vertices[i].y, max.y);
		max.z = fmaxf(mesh->vertices[i].z, max.z);
	}
	bounds = FromMinMax(min, max);
}

mat4 GetWorldMatrix(const Model &model) {

	mat4 translation = Translation(model.position);
	mat4 rotation = Rotation(
			model.rotation.x,
			model.rotation.y,
			model.rotation.z);
	mat4 localMat = rotation * translation;

	mat4 parentMat;
	if (model.parent != 0) {
		parentMat = GetWorldMatrix(*model.parent);
	}
	return localMat * parentMat;
}

OBB GetOBB(const Model &model) {
	mat4 world = GetWorldMatrix(model);
	AABB aabb = model.GetBounds();
	OBB obb;

	obb.size = aabb.size;
	obb.position = MultiplyPoint(aabb.position, world);
	obb.orientation = Cut(world, 3, 3);
	return obb;
}

float ModelRay(const Model &model, const Ray &ray) {
	mat4 world = GetWorldMatrix(model);
	mat4 inv = Inverse(world);

	Ray local;
	local.origin = MultiplyPoint(ray.origin, inv);
	local.direction = MultiplyVector(ray.origin, inv);
	local.NormalizeDirection();

	if (model.GetMesh() != 0) {
		return MeshRay(*(model.GetMesh()), local);
	}
	return -1;
}

bool Linetest(const Model &model, const Line &line) {
	mat4 world = GetWorldMatrix(model);
	mat4 inv = Inverse(world);
	Line local;
	local.start = MultiplyPoint(line.start, inv);
	local.end = MultiplyPoint(line.end, inv);
	if (model.GetMesh() != 0) {
		return Linetest(*(model.GetMesh()), local);
	}
	return false;
}

bool ModelSphere(const Model &model, const Sphere &sphere) {
	mat4 world = GetWorldMatrix(model);
	mat4 inv = Inverse(world);
	Sphere local;
	local.position = MultiplyPoint(sphere.position, inv);
	if (model.GetMesh() != 0) {
		return MeshSphere(*(model.GetMesh()), local);
	}
	return false;
}

bool ModelAABB(const Model &model, const AABB &aabb) {
	mat4 world = GetWorldMatrix(model);
	mat4 inv = Inverse(world);

	OBB local;
	local.size = aabb.size;
	local.position = MultiplyPoint(aabb.position, inv);
	local.orientation = Cut(inv, 3, 3);

	if (model.GetMesh() != 0) {
		return MeshOBB(*(model.GetMesh()), local);
	}
	return false;
}

bool ModelOBB(const Model &model, const OBB &obb) {
	mat4 world = GetWorldMatrix(model);
	mat4 inv = Inverse(world);

	OBB local;
	local.size = obb.size;
	local.position = MultiplyPoint(obb.position, inv);
	local.orientation = obb.orientation * Cut(inv, 3, 3);

	if (model.GetMesh() != 0) {
		return MeshOBB(*(model.GetMesh()), local);
	}
	return false;
}

bool ModelPlane(const Model &model, const Plane &plane) {
	mat4 world = GetWorldMatrix(model);
	mat4 inv = Inverse(world);

	Plane local;
	local.normal = MultiplyVector(plane.normal, inv);
	local.distance = plane.distance;

	if (model.GetMesh() != 0) {
		return MeshPlane(*(model.GetMesh()), local);
	}
	return false;
}

bool ModelTriangle(const Model &model,
		const Triangle &triangle) {

	mat4 world = GetWorldMatrix(model);
	mat4 inv = Inverse(world);

	Triangle local;
	local.a = MultiplyPoint(triangle.a, inv);
	local.b = MultiplyPoint(triangle.b, inv);
	local.c = MultiplyPoint(triangle.c, inv);

	if (model.GetMesh() != 0) {
		return MeshTriangle(*(model.GetMesh()), local);
	}
	return false;
}

Point Intersection(Plane p1, Plane p2, Plane p3) {
	mat3 D{
		p1.normal.x, p2.normal.x, p3.normal.x,
		p1.normal.y, p2.normal.y, p3.normal.y,
		p1.normal.z, p2.normal.z, p3.normal.z
	};

	vec3 A(-p1.distance, -p2.distance, -p3.distance);

	mat3 Dx = D;
	mat3 Dy = D;
	mat3 Dz = D;

	Dx._11 = A.x;
	Dx._12 = A.y;
	Dx._13 = A.z;
	Dy._21 = A.x;
	Dy._22 = A.y;
	Dy._23 = A.z;
	Dz._31 = A.x;
	Dz._32 = A.y;
	Dz._33 = A.z;

	float detD = Determinant(D);

	if (CMP(detD, 0)) {
		return Point();
	}

	float detDx = Determinant(Dx);
	float detDy = Determinant(Dy);
	float detDz = Determinant(Dz);

	return Point(detDx / detD, detDy / detD, detDz / detD);
}

void GetCorners(const Frustum &f, vec3 *outCorners) {
	outCorners[0] = Intersection(f.near, f.top, f.left);
	outCorners[1] = Intersection(f.near, f.top, f.left);
	outCorners[2] = Intersection(f.near, f.bottom, f.left);
	outCorners[3] = Intersection(f.near, f.bottom, f.right);
	outCorners[4] = Intersection(f.far, f.top, f.right);
	outCorners[5] = Intersection(f.far, f.top, f.right);
	outCorners[6] = Intersection(f.far, f.bottom, f.left);
	outCorners[7] = Intersection(f.far, f.bottom, f.right);
}

bool Intersects(const Frustum &f, const Point &p) {
	for (int i = 0; i < 6; ++i) {
		vec3 normal = f.planes[i].normal;
		float dist = f.planes[i].distance;
		float side = Dot(p, normal) + dist;
		if (side < 0.0f) {
			return false;
		}
	}

	return true;
}

bool Intersects(const Frustum &f, const Sphere &s) {
	for (int i = 0; i < 6; ++i) {
		vec3 normal = f.planes[i].normal;
		float dist = f.planes[i].distance;
		float side = Dot(s.position, normal) + dist;

		if (side < -s.radius) {
			return false;
		}
	}

	return true;
}

float Classify(const AABB &aabb, const Plane &plane) {
	float r = fabsf(aabb.size.x * plane.normal.x) + fabsf(aabb.size.y * plane.normal.y) + fabsf(aabb.size.z * plane.normal.z);

	float d = Dot(plane.normal, aabb.position) + plane.distance;

	if (fabsf(d) < r) {
		return 0.0f;
	} else if (d < 0.0f) {
		return d + r;
	}
	return d - r;
}

float Classify(const OBB &obb, const Plane &plane) {
	vec3 normal = MultiplyVector(plane.normal,
			obb.orientation);

	//maximum extent in direction of plane normal
	float r = fabsf(obb.size.x * plane.normal.x) + fabsf(obb.size.y * plane.normal.y) +
			  fabsf(obb.size.z * plane.normal.z);

	float d = Dot(plane.normal, obb.position) + plane.distance;

	if (fabsf(d) < r) {
		return 0.0f;
	} else if (d < 0.0f) {
		return d + r;
	}
	return d - r;
}

bool Intersects(const Frustum &f, const AABB &aabb) {
	for (int i = 0; i < 6; ++i) {
		if (Classify(aabb, f.planes[i]) < 0) {
			return false;
		}
	}
	return true;
}

bool Intersects(const Frustum &f, const OBB &obb) {
	for (int i = 0; i < 6; ++i) {
		if (Classify(obb, f.planes[i]) < 0) {
			return false;
		}
	}
	return true;
}

vec3 Unproject(const vec3 &viewportPoint,
		const vec2 &viewportOrigin, const vec2 &viewportSize,
		const mat4 &view, const mat4 &projection) {
	float normalized[4] = {
		(viewportPoint.x - viewportOrigin.x) / viewportSize.x,
		(viewportPoint.y - viewportOrigin.y) / viewportSize.y,
		viewportPoint.z, 1.0f
	}; //normalized

	float ndcSpace[4] = {
		normalized[0],
		normalized[1],
		normalized[2],
		normalized[3]
	};
	ndcSpace[0] = ndcSpace[0] * 2.0f - 1.0f;

	ndcSpace[1] = 1.0f - ndcSpace[1] * 2.0f;

	if (ndcSpace[2] < 0.0f) {
		ndcSpace[2] = 0.0f;
	}
	if (ndcSpace[2] > 1.0f) {
		ndcSpace[2] = 1.0f;
	}

	mat4 invProjection = Inverse(projection);
	float eyeSpace[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	//eyespace = multiplyPoint(ndcSpace, invProjection);
	Multiply(eyeSpace, ndcSpace, 1, 4,
			invProjection.asArray, 4, 4);

	mat4 invView = Inverse(view);
	float worldSpace[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	//worldspace = multiplypoint(eyespace, invView);
	Multiply(worldSpace, eyeSpace, 1, 4,
			invView.asArray, 4, 4);

	if (!CMP(worldSpace[3], 0.0f)) {
		worldSpace[0] /= worldSpace[3];
		worldSpace[1] /= worldSpace[3];
		worldSpace[2] /= worldSpace[3];
	}

	return vec3(worldSpace[0],
			worldSpace[1],
			worldSpace[2]);
}

Ray GetPickRay(const vec2 &viewportPoint,
		const vec2 &viewportOrigin, const vec2 &viewportSize,
		const mat4 &view, const mat4 &projection) {
	vec3 nearPoint(viewportPoint.x, viewportPoint.y, 0.0f);
	vec3 farPoint(viewportPoint.x, viewportPoint.y, 1.0f);

	vec3 pNear = Unproject(nearPoint, viewportOrigin,
			viewportSize, view, projection);
	vec3 pFar = Unproject(farPoint, viewportOrigin,
			viewportSize, view, projection);
	vec3 normal = Normalized(pFar - pNear);

	vec3 origin = pNear;
	return Ray(origin, normal);
}

void ResetRaycastResult(RaycastResult *outResult) {
	if (outResult != 0) {
		outResult->t = -1;
		outResult->hit = false;
		outResult->normal = vec3(0, 0, 1);
		outResult->point = vec3(0, 0, 0);
	}
}

bool Raycast(const Sphere &sphere,
		const Ray &ray, RaycastResult *outResult) {
	ResetRaycastResult(outResult);
	vec3 e = sphere.position - ray.origin;

	float rSq = sphere.radius * sphere.radius;
	float eSq = MagnitudeSq(e);

	float a = Dot(e, ray.direction);

	float bSq = eSq - (a * a);
	float f = sqrt(rSq - bSq);

	float t = a - f; //Assume normal intersection

	if (rSq - (eSq - a * a) < 0.0f) {
		return false;
	} else if (eSq < rSq) { //Inside sphere
		t = a + f; //Reverse Direction
	}

	if (outResult != 0) {
		outResult->t = t;
		outResult->hit = true;
		outResult->point = ray.origin + ray.direction * t;
		outResult->normal = Normalized(outResult->point -
									   sphere.position);
	}
	return true;
}

void ResetCollisionManifold(CollisionManifold *result) {
	if (result != 0) {
		result->colliding = false;
		result->normal = vec3(0, 0, 1);
		result->depth = FLT_MAX;
		;
		result->contacts.clear();
	}
}

CollisionManifold FindCollisionFeatures(const Sphere &A,
		const Sphere &B) {
	CollisionManifold result;
	ResetCollisionManifold(&result);

	float r = A.radius + B.radius;

	vec3 d = B.position + A.position;
	if (MagnitudeSq(d) - r * r > 0 || MagnitudeSq(d) == 0.0f) {
		return result;
	}
	Normalize(d);

	result.colliding = true;
	result.normal = d;
	result.depth = fabsf(Magnitude(d) - r) * 0.5f;
	//dtp - Distance to intersection point
	float dtp = A.radius - result.depth;
	Point contact = A.position + d * dtp;
	result.contacts.push_back(contact);

	return result;
}

CollisionManifold FindCollisionFeatures(const OBB &A,
		const Sphere &B) {
	CollisionManifold result;
	ResetCollisionManifold(&result);

	Point closestPoint = ClosestPoint(A, B.position);
	float distanceSq = MagnitudeSq(
			closestPoint - B.position);
	if (distanceSq > B.radius * B.radius) {
		return result;
	}

	vec3 normal;
	if (CMP(distanceSq, 0.0f)) {
		float mSq = MagnitudeSq(closestPoint - A.position);
		if (CMP(mSq, 0.0f)) {
			return result;
		}
		//Closest point is at the center of the sphere
		normal = Normalized(closestPoint - A.position);
	} else {
		normal = Normalized(B.position - closestPoint);
	}

	Point outsidePoint = B.position - normal * B.radius;
	float distance = Magnitude(closestPoint - outsidePoint);

	result.colliding = true;
	result.contacts.push_back(closestPoint +
							  (outsidePoint - closestPoint) * 0.05f);

	result.normal = normal;
	result.depth = distance * 0.5f;

	return result;
}

std::vector<Point> GetVertics(const OBB &obb) {
	std::vector<vec3> v;
	v.resize(8);

	vec3 C = obb.position; //OBB Center
	vec3 E = obb.size; //OBB Extents
	const float *o = obb.orientation.asArray;
	vec3 A[] = {
		vec3(o[0], o[1], o[2]),
		vec3(o[3], o[4], o[5]),
		vec3(o[6], o[7], o[8]),
	};

	v[0] = C + A[0] * E[0] + A[1] * E[1] + A[2] * E[2];
	v[1] = C - A[0] * E[0] + A[1] * E[1] + A[2] * E[2];
	v[2] = C + A[0] * E[0] - A[1] * E[1] + A[2] * E[2];
	v[3] = C + A[0] * E[0] + A[1] * E[1] - A[2] * E[2];
	v[4] = C - A[0] * E[0] - A[1] * E[1] - A[2] * E[2];
	v[5] = C + A[0] * E[0] - A[1] * E[1] - A[2] * E[2];
	v[6] = C - A[0] * E[0] + A[1] * E[1] - A[2] * E[2];
	v[7] = C - A[0] * E[0] - A[1] * E[1] + A[2] * E[2];

	return v;
}

std::vector<Line> GetEdges(const OBB &obb) {
	std::vector<Line> result;
	result.reserve(12);

	std::vector<Point> v = GetVertices(obb);

	int index[][2] = { //Indices of edgevertices
		{ 6, 1 }, { 6, 3 }, { 6, 4 }, { 2, 7 }, { 2, 5 }, { 2, 0 },
		{ 0, 1 }, { 0, 3 }, { 7, 1 }, { 7, 4 }, { 4, 5 }, { 5, 3 }
	};

	for (int j = 0; j < 12; ++j) {
		result.push_back(Line(v[index[j][0]], v[index[j][1]]));
	}

	return result;
}

std::vector<Plane> GetPlanes(const OBB &obb) {
	vec3 c = obb.position; //OBBCenter
	vec3 e = obb.size; //OBBExtents
	const float *o = obb.orientation.asArray;

	vec3 a[] = {
		//OBB Axis
		vec3(o[0], o[1], o[2]),
		vec3(o[3], o[4], o[5]),
		vec3(o[6], o[7], o[8]),
	};

	std::vector<Plane> result;
	result.resize(6);

	result[0] = Plane(a[0], Dot(a[0], (c + a[0] * e.x)));
	result[1] = Plane(a[0] * -1.0f, -Dot(a[0], (c - a[0] * e.x)));
	result[2] = Plane(a[1], Dot(a[1], (c + a[1] * e.y)));
	result[3] = Plane(a[1] * -1.0f, -Dot(a[1], (c - a[1] * e.y)));
	result[4] = Plane(a[2], Dot(a[2], (c + a[2] * e.z)));
	result[5] = Plane(a[2] * -1.0f, -Dot(a[2], (c - a[2] * e.z)));

	return result;
}
bool ClipToPlane(const Plane &plane,
		const Line &line, Point *outPoint) {
	vec3 ab = line.end - line.start;
	float nAB = Dot(plane.normal, ab);
	if (CMP(nAB, 0)) {
		return false;
	}

	float nA = Dot(plane.normal, line.start);
	float t = (plane.distance - nA) / nAB;

	if (t >= 0.0f && t < 1.0f) {
		if (outPoint != 0) {
			*outPoint = line.start + ab * t;
		}
		return true;
	}
	return false;
}

std::vector<Point> ClipEdgesToOBB(
		const std::vector<Line> &edges, const OBB &obb) {
	std::vector<Point> result;
	result.reserve(edges.size());
	Point intersection;

	std::vector<Plane> &planes = GetPlanes(obb);

	for (int i = 0; i < planes.size(); ++i) {
		for (int j = 0; j < edges.size(); ++j) {
			if (ClipToPlane(planes[i],
						edges[j], &intersection)) {
				if (PointInOBB(intersection, obb)) {
					result.push_back(intersection);
				}
			}
		}
	}

	return result;
}

float PenetrationDepth(const OBB &o1, const OBB &o2,
		const vec3 &axis, bool *outShouldFlip) {
	Interval i1 = GetInterval(o1, Normalized(axis));
	Interval i2 = GetInterval(o2, Normalized(axis));

	if (!((i2.min <= i1.max) && (i1.min <= i2.max))) {
		return 0.0f; //No penetration
	}

	float len1 = i1.max - i1.min;
	float len2 = i2.max - i2.min;

	float min = fminf(i1.min, i2.min);
	float max = fmaxf(i1.max, i2.max);
	float length = max - min;

	if (outShouldFlip != 0) {
		*outShouldFlip = (i2.min < i1.min);
	}

	return (len1 + len2) - length;
}

CollisionManifold FindCollisionFeatures(const OBB &A,
		const OBB &B) {
	CollisionManifold result;
	ResetCollisionManifold(&result);

	const float *o1 = A.orientation.asArray;
	const float *o2 = B.orientation.asArray;

	vec3 test[15] = { //Face axis
		vec3(o1[0], o1[1], o1[2]),
		vec3(o1[3], o1[4], o1[5]),
		vec3(o1[6], o1[7], o1[8]),
		vec3(o2[0], o2[1], o2[2]),
		vec3(o2[3], o2[4], o2[5]),
		vec3(o2[6], o2[7], o2[8])
	};

	for (int i = 0; i < 3; ++i) { // Fill out rest of axis
		test[6 + i * 3 + 0] = Cross(test[i], test[0]);
		test[6 + i * 3 + 1] = Cross(test[i], test[1]);
		test[6 + i * 3 + 2] = Cross(test[i], test[2]);
	}

	vec3 *hitNormal = 0;
	bool shouldFlip;

	for (int i = 0; i < 15; ++i) {
		if (MagnitudeSq(test[i]) < 0.001f) {
			continue;
		}
		float depth = PenetrationDepth(A, B,
				test[i], &shouldFlip);

		if (depth <= 0.0f) {
			return result;
		} else if (depth < result.depth) {
			if (shouldFlip) {
				test[i] = test[i] * -1.0f;
			}
			result.depth = depth;
			hitNormal = &test[i];
		}
	}

	if (hitNormal == 0) {
		return result;
	}

	vec3 axis = Normalized(*hitNormal);

	std::vector<Point> c1 = ClipEdgesToOBB(GetEdges(B), A);
	std::vector<Point> c2 = ClipEdgesToOBB(GetEdges(A), B);
	result.contacts.reserve(c1.size() + c2.size());
	result.contacts.insert(result.contacts.end(),
			c1.begin(), c1.end());
	result.contacts.insert(result.contacts.end(),
			c2.begin(), c2.end());

	Interval i = GetInterval(A, axis);
	float distance = (i.max - i.min) * 0.5f -
					 result.depth * 0.5f;
	vec3 pointOnPlane = A.position + axis * distance;

	for (int i = result.contacts.size() - 1; i > 0; --i) {
		vec3 contact = result.contacts[i];
		result.contacts[i] = contact + (axis *
											   Dot(axis, pointOnPlane - contact));

		for (int j = result.contacts.size() - 1; j > 1; --j) {
			if (MagnitudeSq(result.contacts[j] - result.contacts[i]) < 0.0001f) {
				result.contacts.erase(result.contacts.begin() + j);
				break;
			}
		}
	}

	result.colliding = true;
	result.normal = axis;

	return result;
}
