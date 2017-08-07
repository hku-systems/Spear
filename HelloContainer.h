#ifndef _HELLOCONTAINER_H
#define _HELLOCONTAINER_H


#include "main.h"

class HelloContainer: public wifiContainer
{
public:
  	int c;
  	int n;

	HelloContainer();
	
	HelloContainer(int id);
	void energySupply1(double theta0);
	void energySupply2();
	void setTheta(double the);	
	void hello_generate_cn();
	void hello_init_output();
	bool isOff(int time);
	bool isOff_collision1(int time);
	void findNext(int time, int timeInterval);
	bool isOff_collision2(int tt, int timeInterval);
};


#endif
