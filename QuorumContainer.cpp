#include "QuorumContainer.h"

using namespace std;
    
QuorumContainer::QuorumContainer(){}

QuorumContainer::QuorumContainer(int id): wifiContainer(id){
    quorum_generate_k();
}

//Energy Module Algorithm 1
void QuorumContainer::energySupply1(double theta0){
	if(energy > Pmin){
		setTheta(theta0 / Pmax * this->energy);
	}
	else{
		setTheta(theta0 / Pmax * Pmin);
	}
}

//Energy Module Algorithm 2
void QuorumContainer::energySupply2(){
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

void QuorumContainer::setTheta(double the){
	theta = the;
	quorum_generate_k();
}
void QuorumContainer::quorum_generate_k(){
	double num;

    for ( int i=1;i<=500;i++ )
    {
	    double val=(2.0*i - 1)/(i*i);
	    if ( val < theta ){
	        double val2=(2*(i-1)-1)/((i-1)*(i-1));
	        if( theta - val < val2 - theta)
	          num=i;
	        else
	          num = i-1;
	        break;
	    }
	}

	for(int i = 0; i < primeset.size(); i++){
		if(primeset[i] == num){
			this->k = primeset[i];
			break;
		}
		else if(primeset[i] > num && i >= 1){
			this->k = primeset[i]-num < num -  primeset[i-1] ? primeset[i] : primeset[i-1];
			break;
		}
		else if(primeset[i] > num && i == 0){
			this->k = primeset[i];
			break;
		}
	}
    rowOn = rand()%k;
    columnOn = rand()%k;
	quorum_init_output();
}

void QuorumContainer::quorum_init_output(){
	stringstream ss;
	ss << "quorum/Prime_Time_pset.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	f << "Node " << ID << ", start at "<<startTime<<", k: " << k <<std::endl;
}

bool QuorumContainer::isOff(int time){
	if( time - startTime > 0){
  		if ( (time - startTime - columnOn) % k == 0 || ( (time - startTime ) % (k*k) >= k*columnOn && (time - startTime ) % (k*k) < k*(columnOn + 1) )){
      		return false;

  		}
	}
	return true;
}
bool QuorumContainer::isOff_collision1(int time){
	int num = rand()%10;
	if(isOff(time) == false){	
		if(num*1.0/10 < PROB1){
			return false;
		}
  	}
	return true;
}

void QuorumContainer::findNext(int time, int timeInterval){
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

bool QuorumContainer::isOff_collision2(int tt, int timeInterval){
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
