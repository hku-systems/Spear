#ifndef _UCONNECTCONTAINER_H
#define _UCONNECTCONTAINER_H

#include "main.h"


class UconnectContainer: public wifiContainer
{
public:
  	int c;
  	int n;
  	
	UconnectContainer();
	UconnectContainer(int id);
	void energySupply1(double theta0);
	void energySupply2();
	void setTheta(double the);	
	void uconnect_generate_cn();
	void uconnect_init_output();
	bool isOff(int time);
	bool isOff_collision1(int time);
	void findNext(int time, int timeInterval);
	bool isOff_collision2(int tt, int timeInterval);
};


#endif
