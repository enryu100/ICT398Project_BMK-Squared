/* register_types.cpp */

#include "register_types.h"
#include "class_db.h"
#include "aiemotion.h"

void register_aiemotion_types(){
	ClassDB::register_class<AIEmotion>();
}

void unregister_aiemotion_types(){
	//do nothing
}