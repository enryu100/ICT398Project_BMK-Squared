/*register_types.cpp */

#include "register_types.h"
#include "class_db.h"
#include "PhysicsSystem.h"

void register_physics_types(){
	ClassDB::register_class<PhysicsSystem>();
}

void unregister_physics_types(){
	
}