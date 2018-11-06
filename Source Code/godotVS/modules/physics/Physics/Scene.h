#pragma once
#ifndef _H_SCENE_
#define _H_SCENE_
#include "Geometry3D.h"
#include <vector>

typedef struct OctreeNode {
	AABB bounds;
	OctreeNode *children;
	std::vector<Model *> models;
	inline OctreeNode() :
		children(0) {}
	inline ~OctreeNode() {
		if (children != 0) {
			delete[] children;
		}
	}
} OctreeNode;

void SplitTree(OctreeNode *node, int depth) {
	if (depth-- <= 0) {//Decrements depth
		return;
	}
	if (node->children == 0) {
		node->children = new OctreeNode[8];
		vec3 c = node->bounds.position;
		vec3 e = node->bounds.size * 0.5f;

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
	if (node->children != 0 && node->models.size() > 0){
		for(int i = 0; i <8; ++i){//for each child
				for (int j = 0, size = node -> models.size();
					j < size; ++j) {
				OBB bounds = GetOBB(*node->models[j]);
					if (AABBOBB(node -> children[i].bounds, bounds)) {
						node->children[i].models.push_back(node->models[j]);
					}
				}
			}
			node->models.clear();

			for (int i = 0; i < 8; ++i) {//recurse
				SplitTree(&(node->children[i]), depth);
			}
		}
}

void Insert(OctreeNode *node, Model *model);
void Remove(OctreeNode *node, Model *model);
void Update(OctreeNode *node, Model *model);

Model *FindClosest(const std::vector<Model *> &set,
		const Ray &ray);
Model *Raycast(OctreeNode *node, const Ray &ray);

std::vector<Model *> Query(OctreeNode *node,
		const Sphere &sphere);
std::vector<Model *> Query(OctreeNode *node,
		const AABB &aabb);

class Scene {
protected:
	std::vector<Model*> objects;
	OctreeNode *octree;

private:
	Scene(const Scene &);
	Scene &operator=(const Scene &);
public:
	inline Scene() :
			octree(0) {}
	inline ~Scene() {
		if (octree != 0) {
			delete octree;
		}
	}

	void AddModel(Model *model);
	void RemoveModel(Model *model);
	void UpdateModel(Model *model);
	std::vector<Model *> FindChildren(const Model *model);

	Model *Raycast(const Ray &ray);
	std::vector<Model *> Query(const Sphere &sphere);
	std::vector<Model *> Query(const AABB &aabb);

	bool Accelerate(const vec3 &position, float size);

	std::vector<Model *> Cull(const Frustum &f);
};

#endif
