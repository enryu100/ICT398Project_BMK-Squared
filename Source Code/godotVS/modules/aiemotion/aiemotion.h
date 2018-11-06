/*AI Emotion Engine Header File*/
#ifndef AIEMOTION_H
#define AIEMOTION_H

#include "reference.h"

class AIEmotion : public Object { 
		GDCLASS(AIEmotion, Object);
	
	protected:
		static void _bind_methods();
	public:
		//increase emotion functions
		int addJoy(double joyVal);
		int addSad(double sadVal);
		int addTrust(double trustVal);
		int addDisgust(double disgVal);
		int addFear(double fearVal);
		int addAnger(double angVal);
		int addSurprise(double surpVal);
		int addAnticipation(double anticVal);
		
		//blank func for registering
		void setNothing();
		
		//getter functions
		double getJoySad(){return joy_sadness;}
		double getTrustDisgust(){return trust_disgust;}
		double getFearAnger(){return fear_anger;}
		double getSurpriseAnticipation(){return surprise_anticipation;}
		
		//setter funcs for modifiable values
		void setEmotionThreshold(double e_threshhold);
		void setDecayRate(r_decay);
		
		//emotion decay function
		void decayEmotions();
	
		//constructor - generic
		AIEmotion();
		// overloaded - sets both in constructor
		AIEmotion(double joy, 
		          double trust,
							double fear,
							double surprise,
							double joy_m, 
							double trust_m,
							double fear_m,
							double surprise_m,
							double sadness_m,
							double disgust_m,
							double anger_m,
							double anticipation_m,
							double e_threshhold);
		//set emotion levels (for use with generic constructor)
		void setEmotionLevels(double joy, 
		          double trust,
							double fear,
							double surprise);
		//set an Agent's modifiers to how quickly they gain points to an emotion
		void setAgentModifiers(double joy_m, 
		                       double trust_m,
                           double fear_m,
                           double surprise_m,
                           double sadness_m,
                           double disgust_m,
                           double anger_m,
                           double anticipation_m);
	private:
		//variable declaration for the emotions
		//default values are 0 (neither one emotion nor the other).
		double joy_sadness,
		       trust_disgust,
           fear_anger,
           surprise_anticipation;
		//variable declaration for the emotion modifiers.
		//linked emotions have separate modifiers as well.
		//defaults are 1, as these are multiplied into the added emotion values
		double m_joy,
           m_sadness,
           m_trust,
           m_disgust,
           m_fear,
           m_anger,
           m_surprise,
           m_anticipation;
		//the emotion level above which an emotional response is triggered
		//default value is 5 (as the scale is from -10 to 10)
		double emotion_threshold;
		//emotion decay rate. multiplied by modifiers.
		double decay_rate;
		//function for normalizing the emotion levels 
		void normalizeEmotions();
}

#endif
