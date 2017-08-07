#ifndef _WIFICONTAINER_H
#define _WIFICONTAINER_H
class wifiContainer
{
public:
	int startTime;
  	long long endTime;
  	double x, y;
  	bool alive;
  	int ID;
  	int energy;
  	double theta;
  	double PROB1;
  	double PROB2;
  	int t1;
  	int t2;
  	bool noOnBefore;
	wifiContainer();
	
	wifiContainer(int id);
	void resetStartTime();
	void reset();
	virtual void energySupply1(double theta0) = 0;
	virtual void energySupply2() = 0;
	virtual void setTheta(double the) = 0;
	virtual bool isOff(int time) = 0;
	virtual	bool isOff_collision1(int time) = 0;
	virtual void findNext(int time, int timeInterval) = 0;
	virtual bool isOff_collision2(int tt, int timeInterval) = 0;
};


#endif
