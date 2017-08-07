#include <string>
#include <sstream>
#include <fstream>
#include "HedisContainer.h"

using namespace std;

HedisContainer::HedisContainer(){}

HedisContainer::HedisContainer(int id): wifiContainer(id){
    hedis_generate_anchor();
}

//Energy Module Algorithm 1
void HedisContainer::energySupply1(double theta0){
	if(energy > Pmin){
		setTheta(theta0 / Pmax * this->energy);
	}
	else{
		setTheta(theta0 / Pmax * Pmin);
	}
}

//Energy Module Algorithm 2
void HedisContainer::energySupply2(){
	for(int i = 0; i < m; i++){
		if(energy <= Pmin){
			setTheta(Pmin*1.0/Pmax * thet[0]);
			break;
		}
		else if(energy > P[i]){
			setTheta(thet[i-1]);
			break;
		}
	}
}

void HedisContainer::setTheta(double the){
	theta = the;
	hedis_generate_anchor();
}

void HedisContainer::hedis_generate_anchor(){
	double num = 2/theta;
	for(int i = 0; i < primeset.size(); i++){
		if(primeset[i] == num){
			this->anchor = primeset[i];
			break;
		}
		else if(primeset[i] > num && i >= 1){
			this->anchor = primeset[i]-num < num -  primeset[i-1] ? primeset[i] : primeset[i-1];
			break;
		}
		else if(primeset[i] > num && i == 0){
			this->anchor = primeset[i];
			break;
		}
	}
	hedis_init_output();
}

void HedisContainer::hedis_init_output(){
	stringstream ss;
	ss << "hedis/Prime_Time_pset.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	f << "Node " << ID << ", start at "<<startTime<<", anchor: " << anchor <<std::endl;
}

bool HedisContainer::isOff(int time){
	if( time - startTime > 0){
  		if ( (time - startTime) % anchor == 0 || (time - startTime - 1) % ( anchor + 1 ) == 0 ){
      		return false;
  		}
	}
	return true;
}

bool HedisContainer::isOff_collision1(int time){
	int num = rand()%10;
	if(isOff(time) == false){
		if(num*1.0/10 < PROB1){
			return false;
		}
  	}
	return true;
}

void HedisContainer::findNext(int time, int timeInterval){
	if(isOff(time) == false){
		t1 = time;
		for(int i = t1+1; i < timeInterval; i++){
			if(isOff(i) == false){
				t2 = i;
				noOnBefore = true;
				return;
			}
		}
		t2 = 0;
	}
}

bool HedisContainer::isOff_collision2(int tt, int timeInterval){
	findNext(tt, timeInterval);
	if(t1 != 0 && t2 != 0 && t1 < t2 && noOnBefore){
		double poss = (t2-tt)*1.0/(t2-t1+1);
		double realP = poss * PROB2;
		int x = rand()%10000;
		if(x <= realP * 10000){
			noOnBefore = false;
			return false;
		}
		return true;
	}
	else if(t1 != 0 && t1 > t2 && t1 == tt){
		return false;
	}
	return true;
}
