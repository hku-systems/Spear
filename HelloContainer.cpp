
#include "HelloContainer.h"

using namespace std;

HelloContainer::HelloContainer(){}

HelloContainer::HelloContainer(int id): wifiContainer(id)
{
    hello_generate_cn();
}

//Energy Module Algorithm 1
void HelloContainer::energySupply1(double theta0){
	if(energy > Pmin){
		setTheta(theta0 / Pmax * this->energy);
	}
	else{
		setTheta(theta0 / Pmax * Pmin);
	}
}

//Energy Module Algorithm 2
void HelloContainer::energySupply2(){
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

void HelloContainer::setTheta(double the){
	theta = the;
	hello_generate_cn();
}

void HelloContainer::hello_generate_cn(){
	double num = 48.5 / ( 97*theta - 1 );
	for(int i = 0; i < primeset.size(); i++){
		if(primeset[i] == num){
			this->n = primeset[i];
			break;
		}
		else if(primeset[i] > num && i >= 1){
			if( primeset[i]-num < num-primeset[i-1]){
				this->n = primeset[i];
			}
			else{
				this->n = primeset[i-1];
			}
			break;
		}
		else if(primeset[i] > num && i == 0){
			this->n = primeset[i];
			break;
		}
	}
	c = 97;
	hello_init_output();
}

void HelloContainer::hello_init_output(){
	std::stringstream ss;
    ss << "hello/Prime Time.txt";
    std::fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
    f << "Node " << ID << ", start at "<<startTime<<", n: "<<n<<", c:"<<c<<std::endl;
}

bool HelloContainer::isOff(int time){
	if(time - startTime > 0){
  		if ( ( time-startTime  ) % (c * n) < c/2.0 || ( time-startTime ) % c == 0){
  			return false;
		}

  		else
    		return true;
	}
	else if (time - startTime == 0)
  		return false;
	else
  		return true;
}
bool HelloContainer::isOff_collision1(int time){
	int num = rand()%10;
	if(isOff(time) == false){	
		if(num*1.0/10 < PROB1){
			return false;
		}
  	}
	return true;
}

void HelloContainer::findNext(int time, int timeInterval){
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

bool HelloContainer::isOff_collision2(int tt, int timeInterval){
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
