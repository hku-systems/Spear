#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>

#include "constant.h"
#include "wifiContainer.h"

wifiContainer::wifiContainer(){}
wifiContainer::wifiContainer(int id)
{
    ID = id;
    alive = true;
    energy = Pmax;
    startTime = rand()%(STARTRANGE+1);
    endTime = -1;
    x = rand()%(SCOPE+1);
    y = rand()%(SCOPE+1);
    theta = 0.1+(rand()%9)*0.05;
	t1 = 0;
    t2 = 0;
    noOnBefore = true;
    PROB1 = 0.5;
    PROB2 = 0.5;
}

void wifiContainer::resetStartTime()
{
	startTime = rand()%(STARTRANGE+1);
}

void wifiContainer::reset()
{
	alive = true;
    energy = Pmax;
}
