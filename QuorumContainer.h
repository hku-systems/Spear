#ifndef _QUORUMCONTAINER_H
#define _QUORUMCONTAINER_H

#include "main.h"

class QuorumContainer: public wifiContainer
{
public:
  	int k; 
    int rowOn;
    int columnOn;
    
	QuorumContainer();
	
	QuorumContainer(int id);
	void energySupply1(double theta0);
	void energySupply2();
	void setTheta(double the);
	void quorum_generate_k();
	void quorum_init_output();
	bool isOff(int time);
	bool isOff_collision1(int time);
	void findNext(int time, int timeInterval);
	bool isOff_collision2(int tt, int timeInterval);
};


#endif
