/*AI Emotion Engine Source File*/

#include "aiemotion.h"
#include <math.h>

void AIEmotion::normalizeEmotions(){
	// check if total means there could be more than 1 value over threshold
  double total;
	total	= fabs(joy_sadness) + fabs(trust_disgust) + fabs(fear_anger) + fabs(surprise_anticipation);
	if(total > (double)19.7){
		//normalization needed
		double ratio = (double)19.7 / total;
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
		return (int)1;
	} else {
		return (int)0;
	}
}

int AIEmotion::addSad(double sadVal){
	joy_sadness -= (sadVal * m_sadness);
	normalizeEmotions();
	if(joy_sadness < (emotion_threshold*-1)){
		return (int)1;
	} else {
		return (int)0;
	}
}

int AIEmotion::addTrust(double trustVal){
	trust_disgust += (trustVal * m_trust);
	normalizeEmotions();
	if(trust_disgust > emotion_threshold){
		return (int)1;
	} else {
		return (int)0;
	}
}

int AIEmotion::addDisgust(double disgVal){
	trust_disgust -= (disgVal * m_disgust);
	normalizeEmotions();
	if(trust_disgust < (emotion_threshold*-1)){
		return (int)1;
	} else {
		return (int)0;
	}
}

int AIEmotion::addFear(double fearVal){
	fear_anger += (fearVal * m_fear);
	normalizeEmotions();
	if(fear_anger > emotion_threshold){
		return (int)1;
	} else {
		return (int)0;
	}
}

int AIEmotion::addAnger(double angVal){
	fear_anger -= (angVal * m_anger);
	normalizeEmotions();
	if(fear_anger < (emotion_threshold*-1)){
		return (int)1;
	} else {
		return (int)0;
	}
}

int AIEmotion::addSurprise(double surpVal){
	surprise_anticipation += (surpVal * m_surprise);
	normalizeEmotions();
	if(surprise_anticipation > emotion_threshold){
		return (int)1;
	} else {
		return (int)0;
	}
}

int AIEmotion::addAnticipation(double anticVal){
	surprise_anticipation -= (anticVal * m_anticipation);
	normalizeEmotions();
	if(surprise_anticipation < (emotion_threshold*-1)){
		return (int)1;
	} else {
		return (int)0;
	}
}

void AIEmotion::calcViewData(){
	int state = 0;
	
	//calculate emotional state based on emotion values
	if (joy_sadness > 0){
		state++;
	} else {
		if(joy_sadness < 0){
			state--;
		}
	}
	if (trust_disgust > 0){
		state++;
	} else {
		if(trust_disgust < 0){
			state--;
		}
	}
	if (fear_anger > 0){
		state++;
	} else {
		if(fear_anger < 0){
			state--;
		}
	}
	if (surprise_anticipation > 0){
		state++;
	} else {
		if(surprise_anticipation < 0){
			state--;
		}
	}
	
	//translate state to view angle and depth
	if(state > 0){
		aov = 120.0;
		dov = 50.0;
	} else if (state < 0){
		aov = 90.0;
		dov = 100.0;
	}
}

void AIEmotion::setEmotionThreshold(double e_threshhold){
	emotion_threshold = e_threshhold;
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

void AIEmotion::setJoyModifier(double joy_m){
	m_joy = joy_m;
}
void AIEmotion::setTrustModifier(double trust_m){
	m_trust = trust_m;
}
void AIEmotion::setFearModifier(double fear_m){
	m_fear = fear_m;
}
void AIEmotion::setSurpriseModifier(double surprise_m){
	m_surprise = surprise_m;
}
void AIEmotion::setSadnessModifier(double sadness_m){
	m_sadness = sadness_m;
}
void AIEmotion::setDisgustModifier(double disgust_m){
	m_disgust = disgust_m;
}
void AIEmotion::setAngerModifier(double anger_m){
	m_anger = anger_m;
}
void AIEmotion::setAnticipationModifier(double anticipation_m){
	m_anticipation = anticipation_m;
}

AIEmotion::AIEmotion(){
	joy_sadness = 0.0;
	trust_disgust = 0.0;
	fear_anger = 0.0;
	surprise_anticipation = 0.0;
	m_joy = 1.0;
	m_sadness = 1.0;
	m_trust = 1.0;
	m_disgust = 1.0;
	m_fear = 1.0;
	m_anger = 1.0;
	m_surprise = 1.0;
	m_anticipation = 1.0;
	emotion_threshold = 5.0;
	decay_rate = 0.05;
	dov = 100.0;
	aov = 120.0;
}

void AIEmotion::_bind_methods(){
	ClassDB::bind_method("addJoy", &AIEmotion::addJoy);
	ClassDB::bind_method("addSad", &AIEmotion::addSad);
	ClassDB::bind_method("addTrust", &AIEmotion::addTrust);
	ClassDB::bind_method("addDisgust", &AIEmotion::addDisgust);
	ClassDB::bind_method("addFear", &AIEmotion::addFear);
	ClassDB::bind_method("addAnger", &AIEmotion::addAnger);
	ClassDB::bind_method("addSurprise", &AIEmotion::addSurprise);
	ClassDB::bind_method("addAnticipation", &AIEmotion::addAnticipation);
	ClassDB::bind_method("setEmotionLevels", &AIEmotion::setEmotionLevels);
	ClassDB::bind_method("setAngerModifier", &AIEmotion::setAngerModifier);
	ClassDB::bind_method("setJoyModifier", &AIEmotion::setJoyModifier);
	ClassDB::bind_method("setFearModifier", &AIEmotion::setFearModifier);
	ClassDB::bind_method("setSadnessModifier", &AIEmotion::setSadnessModifier);
	ClassDB::bind_method("setTrustModifier", &AIEmotion::setTrustModifier);
	ClassDB::bind_method("setDisgustModifier", &AIEmotion::setDisgustModifier);
	ClassDB::bind_method("setSurpriseModifier", &AIEmotion::setSurpriseModifier);
	ClassDB::bind_method("setAnticipationModifier", &AIEmotion::setAnticipationModifier);
	ADD_PROPERTY(PropertyInfo(Variant::REAL, "joy_sadness"), "setNothing","getJoySad");
	ADD_PROPERTY(PropertyInfo(Variant::REAL, "trust_disgust"), "setNothing","getTrustDisgust");
	ADD_PROPERTY(PropertyInfo(Variant::REAL, "fear_anger"), "setNothing","getFearAnger");
	ADD_PROPERTY(PropertyInfo(Variant::REAL, "surprise_anticipation"), "setNothing","getSurpriseAnticipation");
	ADD_PROPERTY(PropertyInfo(Variant::REAL, "e_threshhold"), "setEmotionThreshold", "getEmotionThreshold");
}

void AIEmotion::setNothing(){
	//do nothing
}