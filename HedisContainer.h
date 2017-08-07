#ifndef _HEDISCONTAINER_H
#define _HEDISCONTAINER_H

#include "main.h"


class HedisContainer: public wifiContainer
{
public:
  	int anchor;  //only in hedis
  	
	HedisContainer();
	
	HedisContainer(int id);

	void energySupply1(double theta0);
	void energySupply2();
	void setTheta(double the);
	void hedis_generate_anchor();
	void hedis_init_output();
	bool isOff(int time);
	bool isOff_collision1(int time);
	void findNext(int time, int timeInterval);
	bool isOff_collision2(int tt, int timeInterval);
};

#endif
