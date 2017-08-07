
#include "SearchlightContainer.h"

using namespace std;

SearchlightContainer::SearchlightContainer(){}

SearchlightContainer::SearchlightContainer(int id): wifiContainer(id){
    searchlight_generate_c();
}

//Energy Module Algorithm 1
void SearchlightContainer::energySupply1(double theta0){
	if(energy > Pmin){
		setTheta(theta0 / Pmax * this->energy);
	}
	else{
		setTheta(theta0 / Pmax * Pmin);
	}
}

//Energy Module Algorithm 2
void SearchlightContainer::energySupply2(){
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

void SearchlightContainer::setTheta(double the){
	theta = the;
	searchlight_generate_c();
}
void SearchlightContainer::searchlight_generate_c(){
	double num = 2 / theta;
	for(int i = 0; i < primeset.size(); i++){
		if(primeset[i] == num){
			this->c = primeset[i];
			break;
		}
		else if(primeset[i] > num && i >= 1){
			this->c = primeset[i]-num < num-primeset[i-1] ? primeset[i] : primeset[i-1];
			break;
		}
		else if(primeset[i] > num && i == 0){
			this->c = primeset[i];
			break;
		}
	}
	searchlight_init_output();
}

void SearchlightContainer::searchlight_init_output(){
	std::stringstream ss;
    ss << "searchlight/Prime Time.txt";
    std::fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
    f << "Node " << ID << ", start at "<<startTime<<", c: "<<c<<endl;
}

bool SearchlightContainer::isOff(int time){
	int probe = ceil(c/2);
	int quotient = (time - startTime)/c;
	int remainder = (time - startTime)%c;
	if(time - startTime >= 0){
		if(quotient%probe + 1 == remainder){
			return false;
		}
		if((time - startTime)%c == 0){
			return false;
		}
	}
	return true;
}
bool SearchlightContainer::isOff_collision1(int time){
	int num = rand()%10;
	if(isOff(time) == false){	
		if(num*1.0/10 < PROB1){
			return false;
		}
  	}
	return true;
}

void SearchlightContainer::findNext(int time, int timeInterval){
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

bool SearchlightContainer::isOff_collision2(int tt, int timeInterval){
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
