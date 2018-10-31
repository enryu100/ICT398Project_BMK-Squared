/*AI Emotion Engine Source File*/

#include "aiemotion.h"
#include <math.h>

void AIEmotion::normalizeEmotions(){
	// check if total means there could be more than 1 value over threshold
  double total;
	total	= fabs(joy_sadness) + fabs(trust_disgust) + fabs(fear_anger) + fabs(surprise_anticipation);
	if(total > 19.7d){
		//normalization needed
		double ratio = 19.7d / total;
		joy_sadness *= ratio;
		trust_disgust *= ratio;
		fear_anger *= ratio;
		surprise_anticipation *= ratio;
	}
}

int AIEmotion::addJoy(double joyVal){
	joy_sadness += (joyVal * m_joy);
	normalizeEmotions();
	if(joy_sadness > emotion_threshold){
		return (1);
	} else {
		return (0);
	}
}

int AIEmotion::addSad(double sadVal){
	joy_sadness -= (sadVal * m_sadness);
	normalizeEmotions();
	if(joy_sadness < (emotion_threshold*-1)){
		return (1);
	} else {
		return (0);
	}
}

int AIEmotion::addTrust(double trustVal){
	trust_disgust += (trustVal * m_trust);
	normalizeEmotions();
	if(trust_disgust > emotion_threshold){
		return (1);
	} else {
		return (0);
	}
}

int AIEmotion::addDisgust(double disgVal){
	trust_disgust -= (disgVal * m_disgust);
	normalizeEmotions();
	if(trust_disgust < (emotion_threshold*-1)){
		return (1);
	} else {
		return (0);
	}
}

int AIEmotion::addFear(double fearVal){
	fear_anger += (fearVal * m_fear);
	normalizeEmotions();
	if(fear_anger > emotion_threshold){
		return (1);
	} else {
		return (0);
	}
}

int AIEmotion::addAnger(double angVal){
	fear_anger -= (angVal * m_anger);
	normalizeEmotions();
	if(fear_anger < (emotion_threshold*-1)){
		return (1);
	} else {
		return (0);
	}
}

int AIEmotion::addSurprise(double surpVal){
	surprise_anticipation += (surpVal * m_surprise);
	normalizeEmotions();
	if(surprise_anticipation > emotion_threshold){
		return (1);
	} else {
		return (0);
	}
}

int AIEmotion::addAnticipation(double anticVal){
	surprise_anticipation -= (anticVal * m_anticipation);
	normalizeEmotions();
	if(surprise_anticipation < (emotion_threshold*-1)){
		return (1);
	} else {
		return (0);
	}
}

void AIEmotion::setEmotionThreshold(double e_threshhold){
	emotion_threshold = e_threshhold;
}

void AIEmotion::setDecayRate(r_decay){
	decay_rate = r_decay;
}

void AIEmotion::decayEmotions(){
	//this should be called on a specific interval, to gradually return emotions to a near-midpoint level
	if(joy_sadness > 0){
		joy_sadness *= (decay_rate * m_joy);
	} else {
		joy_sadness *= (decay_rate * m_sadness);
	}
	if(trust_disgust > 0){
		trust_disgust *= (decay_rate * m_trust);
	} else {
		trust_disgust *= (decay_rate * m_disgust);
	}
	if(fear_anger > 0){
		fear_anger *= (decay_rate * m_fear);
	} else {
		fear_anger *= (decay_rate * m_anger);
	}
	if(surprise_anticipation > 0){
		surprise_anticipation *= (decay_rate * m_surprise);
	} else {
		surprise_anticipation *= (decay_rate * m_anticipation);
	}
}

void AIEmotion::setEmotionLevels(double joy, double trust, double fear, double surprise){
	joy_sadness = joy;
	trust_disgust = trust;
	fear_anger = fear;
	surprise_anticipation = surprise;
}

void AIEmotion::setAgentModifiers(double joy_m, double sadness_m, double trust_m, double disgust_m, double fear_m, double anger_m, double surprise_m, double anticipation_m){
	m_joy = joy_m;
	m_sadness = sadness_m;
	m_trust = trust_m;
	m_disgust = disgust_m;
	m_fear = fear_m;
	m_anger = anger_m;
	m_surprise = surprise_m;
	m_anticipation = anticipation_m;
}

AIEmotion::AIEmotion(){
	joy_sadness = 0.0d;
	trust_disgust = 0.0d;
	fear_anger = 0.0d;
	surprise_anticipation = 0.0d;
	m_joy = 1.0d;
	m_sadness = 1.0d;
	m_trust = 1.0d;
	m_disgust = 1.0d;
	m_fear = 1.0d;
	m_anger = 1.0d;
	m_surprise = 1.0d;
	m_anticipation = 1.0d;
	emotion_threshold = 5.0d;
	decay_rate = 0.05d;
}

AIEmotion::AIEmotion(double joy, double trust, double fear, double surprise, double joy_m, double trust_m, double fear_m, double surprise_m, double sadness_m, double disgust_m, double anger_m, double anticipation_m){
	joy_sadness = joy;
	trust_disgust = trust;
	fear_anger = fear;
	surprise_anticipation = surprise;
	m_joy = joy_m;
	m_sadness = sadness_m;
	m_trust = trust_m;
	m_disgust = disgust_m;
	m_fear = fear_m;
	m_anger = anger_m;
	m_surprise = surprise_m;
	m_anticipation = anticipation_m;
	emotion_threshold = 5.0d;
	decay_rate = 0.05d;
}

static void AIEmotion::_bind_methods(){
	ClassDB::register_method(D_METHOD("addJoy", "joyVal"), &addJoy);
	ClassDB::register_method(D_METHOD("addSad", "sadVal"), &addSad);
	ClassDB::register_method(D_METHOD("addTrust", "trustVal"), &addTrust);
	ClassDB::register_method(D_METHOD("addDisgust", "disgVal"), &addDisgust);
	ClassDB::register_method(D_METHOD("addFear", "fearVal"), &addFear);
	ClassDB::register_method(D_METHOD("addAnger", "angVal"), &addAnger);
	ClassDB::register_method(D_METHOD("addSurprise", "surpVal"), &addSurprise);
	ClassDB::register_method(D_METHOD("addAnticipation", "anticVal"), &addAnticipation);
	ClassDB::register_method(D_METHOD("setEmotionLevels", "joy", "trust", "fear", "surprise"), &setEmotionLevels);
	ClassDB::register_method(D_METHOD("setAgentModifiers","joy_m", "sadness_m", "trust_m", "disgust_m", "fear_m", "anger_m", "surprise_m", "anticipation_m"), &setAgentModifiers);
	ClassDB::register_method(D_METHOD("AIEmotion","joy", "trust", "fear", "surprise", "joy_m", "trust_m", "fear_m", "surprise_m", "sadness_m", "disgust_m", "anger_m", "anticipation_m"), &AIEmotion);
	ClassDB::register_method(D_METHOD("AIEmotion"), &AIEmotion);
	ADD_PROPERTY(PropertyInfo(Variant::double, "joy_sadness"), "setNothing","getJoySad");
	ADD_PROPERTY(PropertyInfo(Variant::double, "trust_disgust"), "setNothing","getTrustDisgust");
	ADD_PROPERTY(PropertyInfo(Variant::double, "fear_anger"), "setNothing","getFearAnger");
	ADD_PROPERTY(PropertyInfo(Variant::double, "surprise_anticipation"), "setNothing","getSurpriseAnticipation");
}

void AIEmotion::setNothing(){
	//do nothing
}