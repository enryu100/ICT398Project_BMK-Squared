/*AI Emotion Engine Header File*/
#ifndef AIEMOTION_H
#define AIEMOTION_H

#include "reference.h"

/********************************************//**
 * \class AIEmotion aiemotion.h "aiemotion.h"
 * \brief Manages the AI's emotion system.
 * \details Contains the variables required to simulate an NPC's emotions, based on Plutchik's Wheel of Emotions, and the functions required to manipulate them
 * \author Brandon Jin Yang Lim
 * \version 1.0.10a
 * \date 30-10-2018
 * \pre Should be declared as a variable within the NPC class.
 * \bug Does not support simulation of multiple high emotions at once
 * \todo Ensure the emotional decay function works and can be implemented
 ***********************************************/
class AIEmotion : public Object { 
		GDCLASS(AIEmotion, Object);
	
	protected:
		static void _bind_methods();
	public:
		//increase emotion functions
		/********************************************//**
		 * \defgroup emotionIncrease Emotion value increase functions
		 * @{
		 ***********************************************/
		/********************************************//**
		 * \fn addJoy
		 * \brief Adds to the value of Joy
		 * \details Parameter passed is added to the value of joy (reducing sadness), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \param[in] joyVal Amount to increase joy by.
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addJoy(double joyVal);
		
		/********************************************//**
		 * \fn addJoy
		 * \brief Adds to the value of Sadness
		 * \details Parameter passed is added to the value of sadness (reducing joy), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \param[in] sadVal Amount to increase sadness by.
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addSad(double sadVal);
		
		/********************************************//**
		 * \fn addTrust
		 * \brief Adds to the value of Trust
		 * \details Parameter passed is added to the value of trust (reducing disgust), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \param[in] trustVal Amount to increase trust by.
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addTrust(double trustVal);
		
		/********************************************//**
		 * \fn addDisgust
		 * \brief Adds to the value of Disgust
		 * \details Parameter passed is added to the value of disgust (reducing trust), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \param[in] disgVal Amount to increase disgust by.
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addDisgust(double disgVal);
		
		/********************************************//**
		 * \fn addFear
		 * \brief Adds to the value of Fear
		 * \details Parameter passed is added to the value of fear (reducing anger), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \param[in] fearVal Amount to increase fear by.
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addFear(double fearVal);
		
		/********************************************//**
		 * \fn addAnger
		 * \brief Adds to the value of Anger
		 * \details Parameter passed is added to the value of anger (reducing fear), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addAnger(double angVal);
		
		/********************************************//**
		 * \fn addSurprise
		 * \brief Adds to the value of Surprise
		 * \details Parameter passed is added to the value of surprise (reducing anticipation), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \param[in] surpVal Amount to increase surprise by.
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addSurprise(double surpVal);
		
		/********************************************//**
		 * \fn addAnticipation
		 * \brief Adds to the value of Anticipation
		 * \details Parameter passed is added to the value of anticipation (reducing surprise), then emotions are normalized. Returns a 1 if the emotion value passes a given threshold
		 * \param[in] anticVal Amount to increase anticipation by.
		 * \author Brandon Jin Yang Lim
		 * \return flag showing whether emotion has passed threshold
		 ***********************************************/
		int addAnticipation(double anticVal);
		/**
		 * @} //end of EmotionIncrease group
		 */
		
		//blank func for registering
		/********************************************//**
		 * \fn setNothing
		 * \brief Placeholder function for Godot module config.
		 * \details Empty function to act as a placeholder for declaring variables in Godot. Godot requires getters and setters for any variable declared.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setNothing();
		
		//getter functions
		/********************************************//**
		 * \defgroup getters Get Value functions
		 * @{
		 ***********************************************/
		/********************************************//**
		 * \fn getJoySad
		 * \brief gets the value on the Joy-Sadness axis
		 * \author Brandon Jin Yang Lim
		 * \return value of joy_sadness
		 ***********************************************/
		double getJoySad(){return joy_sadness;}
		
		/********************************************//**
		 * \fn getTrustDisgust
		 * \brief gets the value on the Trust-Disgust axis
		 * \author Brandon Jin Yang Lim
		 * \return value of trust_disgust
		 ***********************************************/
		double getTrustDisgust(){return trust_disgust;}
		
		/********************************************//**
		 * \fn getFearAnger
		 * \brief gets the value on the Fear-Anger axis
		 * \author Brandon Jin Yang Lim
		 * \return value of fear_anger
		 ***********************************************/
		double getFearAnger(){return fear_anger;}
		
		/********************************************//**
		 * \fn getSurpriseAnticipation
		 * \brief gets the value on the Surprise-Anticipation axis
		 * \author Brandon Jin Yang Lim
		 * \return value of surprise_anticipation
		 ***********************************************/
		double getSurpriseAnticipation(){return surprise_anticipation;}
		
		/********************************************//**
		 * \fn getEmotionThreshold
		 * \brief gets the emotionThreshold value
		 * \author Brandon Jin Yang Lim
		 * \return value of emotion_threshold
		 ***********************************************/
		double getEmotionThreshold(){return emotion_threshold;}
		
		/********************************************//**
		 * \fn getDepthOfView
		 * \brief gets the depth of view of the NPC
		 * \details depth of view is calculated based on emotional state
		 * \author Brandon Jin Yang Lim
		 * \return depth of view, in units
		 ***********************************************/
		double getDepthOfView(){return dov;}
		
		/********************************************//**
		 * \fn getAngleOfView
		 * \brief gets the angle of view of the NPC
		 * \author Brandon Jin Yang Lim
		 * \return angle of view, in degrees
		 ***********************************************/
		 double getAngleOfView(){return aov;}
		/********************************************//**
		 * @}
		 ***********************************************/
		 
		//calculate view angle/depth
		/********************************************//**
		 * \fn calcViewData
		 * \brief calculates the depth and angle of view
		 * \details calculates depth and angle of view based on emotional state.
		 * \todo use actual calculations instead of two sets of arbitrary values.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void calcViewData();
		
		//emotion decay function
		/********************************************//**
		 * \fn decayEmotions
		 * \brief reduces emotions towards their axes' centres
		 * \details divides emotion values by a set amount whenever called.
		 * \pre should be called within a timer function
		 * \warning beta version. does not work perfectly.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void decayEmotions();
		
		/********************************************//**
		 * \defgroup setters Setter functions
		 * @{
		 ***********************************************/
		
		//setter funcs for modifiable values
		/********************************************//**
		 * \fn setEmotionThreshold
		 * \brief sets threshold at which response is triggered
		 * \param[in] emotion_threshold Value to which the emotion threshold is to be set.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setEmotionThreshold(double emotion_threshold);
		
		//set emotion levels (for use with generic constructor)
		/********************************************//**
		 * \fn setEmotionLevels
		 * \brief sets the values for each emotion axis
		 * \param[in] joy Intended value of joy_sadness.
		 * \param[in] trust Intended value of trust_disgust.
		 * \param[in] fear Intended value of fear_anger.
		 * \param[in] surprise Intended value of surprise_anticipation.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setEmotionLevels(double joy, 
		          double trust,
							double fear,
							double surprise);
		
		//set an Agent's modifiers to how quickly they gain points to an emotion
		/********************************************//**
		 * \fn setJoyModifier
		 * \brief sets the modifier for joy input
		 * \param[in] joy_m Intended joy modifier.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setJoyModifier(double joy_m);
		
		/********************************************//**
		 * \fn setTrustModifier
		 * \brief sets the modifier for trust input
		 * \param[in] trust_m Intended trust modifier.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setTrustModifier(double trust_m);
		
		/********************************************//**
		 * \fn setFearModifier
		 * \brief sets the modifier for fear input
		 * \param[in] fear_m Intended fear modifier.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setFearModifier(double fear_m);
		
		/********************************************//**
		 * \fn setSurpriseModifier
		 * \brief sets the modifier for surprise input
		 * \param[in] surprise_m Intended surprise modifier.
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setSurpriseModifier(double surprise_m);
		
		/********************************************//**
		 * \fn setSadnessModifier
		 * \brief sets the modifier for sadness input
		 * \param[in] sadness_m Intended sadness modifier.
		 * \note Is separate to joy modifier to be more realistic
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setSadnessModifier(double sadness_m);
		
		/********************************************//**
		 * \fn setDisgustModifier
		 * \brief sets the modifier for disgust input
		 * \param[in] disgust_m Intended disgust modifier.
		 * \note Is separate to trust modifier to be more realistic
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setDisgustModifier(double disgust_m);
		
		/********************************************//**
		 * \fn setAngerModifier
		 * \brief sets the modifier for anger input
		 * \param[in] anger_m Intended anger modifier.
		 * \note Is separate to fear modifier to be more realistic
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setAngerModifier(double anger_m);
		
		/********************************************//**
		 * \fn setAnticipationModifier
		 * \brief sets the modifier for anticipation input
		 * \param[in] anticipation_m Intended anticipation modifier.
		 * \note Is separate to surprise modifier to be more realistic
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void setAnticipationModifier(double anticipation_m);
		
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
							double emotion_threshold);
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
					 
		// angle/depth of view variables. depth in "units", angle in degrees
		double dov, aov;
		//the emotion level above which an emotional response is triggered
		//default value is 5 (as the scale is from -10 to 10)
		double emotion_threshold;
		//emotion decay rate. multiplied by modifiers.
		double decay_rate;
		//function for normalizing the emotion levels 
		/********************************************//**
		 * \fn normalizeEmotions
		 * \brief ensures emotion values don't exceed maximums
		 * \details Manipulates the emotion variables to ensure only 1 value reaches the maximum at a time, and no value can be excessively high
		 * \author Brandon Jin Yang Lim
		 ***********************************************/
		void normalizeEmotions();
};

#endif
