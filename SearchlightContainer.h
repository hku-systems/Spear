#ifndef _SEARCHLIGHTCONTAINER_H
#define _SEARCHLIGHTCONTAINER_H

#include "main.h"


class SearchlightContainer: public wifiContainer
{
public:
  	int c;
  	
	SearchlightContainer();
	
	SearchlightContainer(int id);
	
	void energySupply1(double theta0);
	void energySupply2();
	void setTheta(double the);
	void searchlight_generate_c();
	void searchlight_init_output();
	bool isOff(int time);
	bool isOff_collision1(int time);
	void findNext(int time, int timeInterval);
	bool isOff_collision2(int tt, int timeInterval);
};

#endif
