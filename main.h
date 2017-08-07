#ifndef _MAIN_H
#define _MAIN_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>


#include "wifiContainer.h"
#include "constant.h"

extern int latencyMatrix[1000][1000];
extern int neighborMap[1000][1000];
extern int aveLatency[1000][1000];
extern std::vector<int> primeset;
extern int P[m];
extern double thet[m];
extern double dutycycle[TIMES];
extern double percent[TIMES];
extern int late[TIMES];
extern wifiContainer *warray[1000];
extern wifiContainer *w[NODENUMBOUND];
extern int latePUT[NODENUMPUT][NODENUMPUT];
extern double neighborMapPUT[NODENUMPUT][NODENUMPUT];


#endif
