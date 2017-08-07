/* 
 * Created by Zhen Cao on Aug 6 2017
 * Codes for simulation of Spear: a practical neighbor discovery framework for wireless sensor networks
 * If there is any problem, please contact caozhen_bupt@hotmail.com
 * Due to the limitations of time, I didn't consider the invalid input. Please input as instructed.
 * It can be used in non-commercial situation. 
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>


#include "HedisContainer.h"
#include "HelloContainer.h"
#include "SearchlightContainer.h"
#include "UconnectContainer.h"
#include "QuorumContainer.h"

/* infinite number
 */
#define INF    0x3f3f3f3f

using namespace std;

/* Matrix for 1000 nodes
 */
int latencyMatrix[1000][1000];
int neighborMap[1000][1000];
int aveLatency[1000][1000];

/* a vector to store primes
 */
std::vector<int> primeset;

/* PWR's P and theta sequence
 */
int P[m];
double thet[m];
/* sensor nodes for 1000 nodes and star network
 */
wifiContainer *warray[1000];
wifiContainer *w[NODENUMBOUND];

/* array to store dutycycle, energy percentage, 
 * latency in energy part: Bare Protocao, PWR, PDR
 */
double dutycycle[TIMES];
double percent[TIMES];
int late[TIMES];

/*Matrix for throughput: latency and neighbor map*/
int latePUT[NODENUMPUT][NODENUMPUT];
double neighborMapPUT[NODENUMPUT][NODENUMPUT];


/* input: number n
 * output: vector<int> primeset 
 * function: calculate all primes that smaller than n
 */
void findAllPrimes(int n)
{
    primeset.push_back(2);
    for(int k = 3; k<=n; k++){
        bool isPrime = true;
        for(int i = 2; i*i <= k; i++){
            if(k % i == 0){
                isPrime = false;
                break;
            }
        }
        if(isPrime){
            primeset.push_back(k);
        }
    }
}


/* input: theta0: original dutycycle for PDR; dir: directory to store result files
 * output: P[i], thet[i]
 * function: calculate remaining energy and dutycycle bound in advance for star networks
 */
void initEnergySupply2Star(double theta0, string dir){
	stringstream ss;
	ss <<dir<< "/star_energy_theta_table.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	w[0]->energy = Pmax;
	w[0]->setTheta(theta0);
	w[0]->alive = true;
	P[0] = Pmax;
	thet[0] = theta0;
	for(int i = 1; i < m; i++){
		P[i] = P[i-1] - T * thet[i-1] * pn;
		thet[i] = P[i]*1.0/P[0] * thet[0];
		f<<i<<" "<<P[i]<<" "<<thet[i]<<endl;
	}
}


/* input: theta0: original dutycycle for PDR; dir: directory to store result files
 * output: P[i], thet[i]
 * function: calculate remaining energy and dutycycle bound in advance for 100*100 throughput network
 */
void initEnergySupply2PUT(double theta0, string dir){
	stringstream ss;
	ss <<dir<< "/throughput_energy_theta_table.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	P[0] = Pmax;
	thet[0] = theta0;
	for(int i = 1; i < m; i++){
		P[i] = P[i-1] - T * thet[i-1] * pn;
		thet[i] = P[i]*1.0/P[0] * thet[0];
		f<<i<<" "<<P[i]<<" "<<thet[i]<<endl;
	}
}


/* initialize throughput network latencymatrix latePUT
 */
void initLatePUT()
{
	for(int i = 0; i < NODENUMPUT; i++){
		for (int j = 0; j < NODENUMPUT; j++){
      		latePUT[i][j] = INF;
    	}
  	}
}


/* set throughput network latencymatrix latePUT,
 * node i and j discover each other at time slot t
 */
void setLatePUT(int i, int j, int t)
{
  	if(latePUT[i][j] == INF)
  	{
    	if (warray[i]->startTime > warray[j]->startTime){
    		latePUT[i][j] = t - warray[i]->startTime;
    		latePUT[j][i] = t - warray[i]->startTime;
    		
    	}

    	else{
    		latePUT[i][j] = t - warray[j]->startTime;
    		latePUT[j][i] = t - warray[j]->startTime;
    	}
  	}
}


/* initialize throughput network nodes' adjacent matrix neighborMapPUT,
 * if node i and j are neighbors to each other, neighborMapPUT[i][j]=1
 * else neighborMapPUT[i][j]=INF
 */
string initMapPUT(int alg)
{
	string dir;
	switch(alg){
		case 1:{
			dir = "hedis";
			for(int i = 0; i < NODENUMPUT; i++){
				warray[i] =  new HedisContainer(i);
			}
			break;
		}
		case 2:{
			dir = "hello";
			for(int i = 0; i < NODENUMPUT; i++){
				warray[i] =  new HelloContainer(i);
			}
			break;
		}
		case 3:{
			dir = "searchlight";
			for(int i = 0; i < NODENUMPUT; i++){
				warray[i] =  new SearchlightContainer(i);
			}
			break;
		}
		case 4:{
			dir = "uconnect";
			for(int i = 0; i < NODENUMPUT; i++){
				warray[i] =  new UconnectContainer(i);
			}
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			dir="";
			return dir;
	}
	stringstream ss;
	ss <<dir<< "/throughput_neighborMap.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	for(int i = 0; i <  NODENUMPUT; i++){
		for(int j = 0; j <NODENUM; j++){
			if(i==j){
				neighborMapPUT[i][j] = 0;
			}
			else{
				neighborMapPUT[i][j] = INF;
			}
		}
	}
	for(int i = 0; i < NODENUMPUT; i++){
		for(int j = 0; j < NODENUMPUT; j++){
			double tmp = sqrt((warray[i]->x-warray[j]->x)*(warray[i]->x-warray[j]->x)+(warray[i]->y-warray[j]->y)*(warray[i]->y-warray[j]->y));
			if(tmp <= dc && i != j){
				neighborMapPUT[i][j] = 1;
				neighborMapPUT[j][i] = 1;
			}
		}
	}
	for(int i = 0; i < NODENUMPUT; i++){
		for(int j = 0; j < NODENUMPUT; j++){
			f<<neighborMapPUT[i][j]<<",";
		}
		f<<endl;
	}
	return dir;
}


/* initialize the latency matrix of NODENUMPUT nodes 
 * distribution network: latencyMatrix
 */
void initLatency()
{
	for(int i = 0; i < 1000; i++){
		for (int j = 0; j < 1000; j++){
      		latencyMatrix[i][j] = MAX_INT;
    	}
  	}
}



/* set 1000 nodes distribution network latencymatrix latencyMatrix,
 * node i and j discover each other at time slot t
 */
void setLatency(int i, int j, int t)
{
  	if(latencyMatrix[i][j] == MAX_INT)
  	{
    	if (warray[i]->startTime > warray[j]->startTime){
    		latencyMatrix[i][j] = t - warray[i]->startTime;
    		latencyMatrix[j][i] = t - warray[i]->startTime;
    		
    	}

    	else{
    		latencyMatrix[i][j] = t - warray[j]->startTime;
    		latencyMatrix[j][i] = t - warray[j]->startTime;
    	}
  	}
}



/* initialize the neighbor map in throughput
 */
string initMapDistr(int alg)
{
	string dir;
	switch(alg){
		case 1:{
			dir = "hedis";
			for(int i = 0; i < 1000; i++){
				warray[i] =  new HedisContainer(i);
			}
			break;
		}
		case 2:{
			dir = "hello";
			for(int i = 0; i < 1000; i++){
				warray[i] =  new HelloContainer(i);
			}
			break;
		}
		case 3:{
			dir = "searchlight";
			for(int i = 0; i < 1000; i++){
				warray[i] =  new SearchlightContainer(i);
			}
			break;
		}
		case 4:{
			dir = "uconnect";
			for(int i = 0; i < 1000; i++){
				warray[i] =  new UconnectContainer(i);
			}
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			dir="";
			return dir;
	}
	

	for(int i = 0; i < 1000; i++){
		for(int j = 0; j < 1000; j++){
			if(sqrt((warray[i]->x-warray[j]->x)*(warray[i]->x-warray[j]->x)+(warray[i]->y-warray[j]->y)*(warray[i]->y-warray[j]->y)) <= dc && i != j){//warray[i]->alive == true && warray[j]->alive == true &&
				neighborMap[i][j] = 1;
				neighborMap[j][i] = 1;
			}
		}
	}
	stringstream nb;
	nb <<dir<< "/distribution_neighbormap.txt";
	fstream nbf(nb.str().c_str(), std::ios::out | std::ios::app);
	for(int i = 0; i < 1000; i++){
		for(int j = 0; j < 1000; j++){
			nbf<<neighborMap[i][j]<<",";
		}
		nbf<<endl;
	}
}


/* Input: timeInterval: the number of time slots in a time period
		  arrayLength: to calculate arrayLength times average number
		  alg: which algorithm
	Output: duty cycle, average latency, max latency
	Function: calculate average latency and max latency for two nodes in 
	symmetric duty cycle
 */
void twoNodeSym(int timeInterval, int arrayLength, int alg)
{
	int twoNodeAvgLatency[21]={0};
	int twoNodeMaxLatency[21]={0};
	string dir;
	wifiContainer *A, *B;
	switch(alg){
		case 1:{
			dir = "hedis";
			A = new HedisContainer(0);
			B = new HedisContainer(1);
			break;
		}
		case 2:{
			dir = "hello";
			A = new HelloContainer(0);
			B = new HelloContainer(1);
			break;
		}
		case 3:{
			dir = "searchlight";
			A = new SearchlightContainer(0);
			B = new SearchlightContainer(1);
			break;
		}
		case 4:{
			dir = "uconnect";
			A = new UconnectContainer(0);
			B = new UconnectContainer(1);
			break;
		}
		case 5:{
			dir = "quorum";
			A = new QuorumContainer(0);
			B = new QuorumContainer(1);
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
	
	for (int i = 0; i < arrayLength; i++){
	  	int latency;
	  	double timeFrame;

	  	A->startTime = rand()%(STARTRANGE+1);
	  	B->startTime = rand()%(STARTRANGE+1);
	  	int laterTime = A->startTime > B->startTime ? A->startTime : B->startTime;
	  	for(double bthe = 0.1; bthe <= 0.31; bthe += 0.01) {
	  		latency = MAX_INT;
	      	A->setTheta(bthe);
	  		B->setTheta(bthe);
	      	A->resetStartTime();
	  		B->resetStartTime();
	  		laterTime = A->startTime > B->startTime ? A->startTime : B->startTime;
	
	  		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame++){
	  			if(A->isOff(timeFrame) == false && B->isOff(timeFrame) == false){//i and j are neighbors in default
	  				latency = timeFrame - laterTime;
	  				break;
	  			}
	  		}
	    	double temp = (bthe - 0.10)/0.01;
	      	int nodeNum = (int) (temp+0.5);
	      	if( latency != MAX_INT){
	        	twoNodeAvgLatency[nodeNum] += latency;
	        	if (twoNodeMaxLatency[nodeNum] < latency)
	          		twoNodeMaxLatency[nodeNum] = latency;
	      	}
	    }
  	}
    stringstream ss;
    ss << dir<<"/two_node_latency_sym.txt";
    fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
    for (int i=0; i<21; i++) {
	    double th = i*0.01 + 0.1;
	    twoNodeAvgLatency[i]=twoNodeAvgLatency[i]*1.0/arrayLength;
	    f <<th << "," << twoNodeAvgLatency[i] <<","<< twoNodeMaxLatency[i] <<endl;
    }
}



/* Input: timeInterval: the number of time slots in a time period
		  arrayLength: to calculate arrayLength times average number
		  alg: which algorithm
	Output: duty cycle, average latency, max latency
	Function: calculate average latency and max latency for two nodes in 
	asymmetric duty cycle
 */
void twoNodeAsym(int timeInterval, int arrayLength, int alg)
{
	int twoNodeAvgLatency[21]={0};
	int twoNodeMaxLatency[21]={0};
	string dir;
	wifiContainer *A, *B;
	switch(alg){
		case 1:{
			dir = "hedis";
			A = new HedisContainer(0);
			B = new HedisContainer(1);
			break;
		}
		case 2:{
			dir = "hello";
			A = new HelloContainer(0);
			B = new HelloContainer(1);
			break;
		}
		case 3:{
			dir = "searchlight";
			A = new SearchlightContainer(0);
			B = new SearchlightContainer(1);
			break;
		}
		case 4:{
			dir = "uconnect";
			A = new UconnectContainer(0);
			B = new UconnectContainer(1);
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
	
	for (int i = 0; i < arrayLength; i++) {
	  	A->setTheta(0.2);
	  	int latency;
	  	double timeFrame;
	  	A->startTime = rand()%(STARTRANGE+1);
	  	B->startTime = rand()%(STARTRANGE+1);
	  	int laterTime = A->startTime > B->startTime ? A->startTime : B->startTime;

	  	for(double bthe = 0.1; bthe <= 0.31; bthe += 0.01) {
	  		latency = MAX_INT;
	  		B->setTheta(bthe);
	  		B->resetStartTime();
	  		laterTime = A->startTime > B->startTime ? A->startTime : B->startTime;
	
	  		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame ++){
	  			if(A->isOff(timeFrame) == false && B->isOff(timeFrame) == false){//i and j are neighbors in default
	  				latency = timeFrame - laterTime;
	  				break;
	  			}
	  		}
	      	double temp = (bthe - 0.10)/0.01;
	      	int nodeNum = (int) (temp+0.5);
	      	if( latency != MAX_INT){
		        twoNodeAvgLatency[nodeNum] += latency;
		        if (twoNodeMaxLatency[nodeNum] < latency)
		          twoNodeMaxLatency[nodeNum] = latency;
	    	}
    	}
  	}
    stringstream ss;
    ss << dir<<"/two_node_latency_asym.txt";
    fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
    for (int i=0; i<21; i++) {
	    double th = i*0.01 + 0.1;
	    twoNodeAvgLatency[i]=twoNodeAvgLatency[i]/arrayLength*1.0;
	    f <<th << "," << twoNodeAvgLatency[i] <<","<< twoNodeMaxLatency[i] <<endl;
    }
}



/*  Input: timeInterval: the number of time slots in a time period
 *		  alg: which algorithm
 *	Output: duty cycle, remaining energy percentage, average latency after each timeInterval
 *	Function: energy percentage and life time of a bare protocol in star networks
 */
void energy_0(int timeInterval, int alg)
{
	int findNeighborTime[NODENUM];
	int onNeighbor = 0;
	int max = 0;
	double theta0 = 0.3;
	int flag = true;
	int bound = 0;
	string dir = "";
	
	for(int i = 0; i < NODENUM; i++){
		findNeighborTime[i] = 0;
	}
	
	switch(alg){
		case 1:{
			dir = "hedis";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HedisContainer(i);
			}
			break;
		}
		case 2:{
			dir = "hello";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HelloContainer(i);
			}			
			break;
		}
		case 3:{
			dir = "searchlight";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new SearchlightContainer(i);
			}
			break;
		}
		case 4:{
			dir = "uconnect";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new UconnectContainer(i);
			}
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
	
	stringstream ss;
	ss << dir<<"/star_energy_notchange.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	w[0]->setTheta(theta0);
	w[0]->resetStartTime();
	w[0]->endTime = 0;
	for(int k = 0; k < TIMES && w[0]->alive == true; k++){
		bound = k;
		for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval && w[0]->alive == true; timeFrame ++){
			w[0]->endTime ++;
			if(w[0]->isOff(timeFrame) == false  && w[0]->alive == true){
				w[0]->energy -= pn;
				if(w[0]->energy <= 0 && flag == true){
					w[0]->alive = false;
					flag = false;
					break;
				}
				for(int j = 1; j < NODENUM; j++){
					if(w[j]->isOff(timeFrame) == false && w[j]->alive == true){
						onNeighbor++;
					}
				}
				if(onNeighbor == 1){
					for(int j = 1; j < NODENUM; j++){
						if(w[j]->isOff(timeFrame) == false && findNeighborTime[j] == 0 && w[j]->alive == true){
							findNeighborTime[j] = timeFrame;
							max = timeFrame;
						}
					}
				}
			}
			onNeighbor = 0;
		}
		
		dutycycle[k] += w[0]->theta;
		percent[k] += w[0]->energy*1.0 / Pmax;
		
		for(int l = 1; l < NODENUM; l++){
			findNeighborTime[l] = 0;
		}
		if(max == 0){
			late[k] += timeInterval - w[0]->startTime;
		}
		else{
			late[k] += max - w[0]->startTime;
		}
		max = 0;
	}
	
	for(int k = 0; k < bound; k++){
		f<<dutycycle[k]<<","<<percent[k]<<","<<late[k]<<endl;
		
	}
	f<<"life Cycle = "<<w[0]->endTime - w[0]->startTime<<endl;
}



/*  Input: timeInterval: the number of time slots in a time period
 *		  alg: which algorithm
 *	Output: duty cycle, remaining energy percentage, average latency after each timeInterval
 *	Function: energy percentage and life time of PDR in star networks
 */
void energy_1(int timeInterval, int alg)
{	
	int findNeighborTime[NODENUM];
	int onNeighbor = 0;
	int max = 0;
	double theta0 = 0.3;
	int flag = true;
	string dir = "";
	int bound = 0;
	
	for(int i = 0; i < NODENUM; i++){
		findNeighborTime[i] = 0;
	}
	switch(alg){
		case 1:{
			dir = "hedis";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HedisContainer(i);
			}
			break;
		}
		case 2:{
			dir = "hello";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HelloContainer(i);
			}			
			break;
		}
		case 3:{
			dir = "searchlight";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new SearchlightContainer(i);
			}
			break;
		}
		case 4:{
			dir = "uconnect";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new UconnectContainer(i);
			}
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
		
	stringstream ss;
	ss <<dir<< "/star_energy1.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	w[0]->setTheta(theta0);
	w[0]->resetStartTime();
	for(int k = 0; k < TIMES && w[0]->alive == true; k++){
		bound = k;
		for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval && w[0]->alive == true; timeFrame ++){
			w[0]->endTime ++;
			if(w[0]->isOff(timeFrame) == false  && w[0]->alive == true){
				w[0]->energy -= pn;
				if(w[0]->energy <= 0 && flag == true){
					w[0]->alive = false;
					flag = false;
					break;
				}
				for(int j = 1; j < NODENUM; j++){
					if(w[j]->isOff(timeFrame) == false && w[j]->alive == true){
						onNeighbor++;
					}
				}
				if(onNeighbor == 1){
					for(int j = 1; j < NODENUM; j++){
						if(w[j]->isOff(timeFrame) == false && findNeighborTime[j] == 0 && w[j]->alive == true){
							findNeighborTime[j] = timeFrame;
							max = timeFrame;
						}
					}
				}
			}
			onNeighbor = 0;
		}
		
		w[0]->energySupply1(theta0);
		dutycycle[k] += w[0]->theta;
		percent[k] += w[0]->energy*1.0 / Pmax;
		for(int l = 1; l < NODENUM; l++){
			findNeighborTime[l] = 0;
		}
		if(max == 0){
			late[k] += timeInterval - w[0]->startTime;
		}
		else{
			late[k] += max - w[0]->startTime;
		}
		max = 0;
	}
	for(int k = 0; k < bound; k++){
		f<<dutycycle[k]<<","<<percent[k]<<","<<late[k]<<endl;
		
	}
	f<<"life Cycle = "<<w[0]->endTime - w[0]->startTime<<endl;
}



/*  Input: timeInterval: the number of time slots in a time period
 *		  alg: which algorithm
 *	Output: duty cycle, remaining energy percentage, average latency after each timeInterval
 *	Function: energy percentage and life time of PWR in star networks
 */
void energy_2(int timeInterval, int alg)
{
	int findNeighborTime[NODENUM];
	int onNeighbor = 0;
	int max = 0;
	double theta0 = 0.3;
	int flag = true;
	string dir = "";
	int bound = 0;
	
	for(int i = 0; i < NODENUM; i++){
		findNeighborTime[i] = 0;
	}
	switch(alg){
		case 1:{
			dir = "hedis";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HedisContainer(i);
			}
			break;
		}
		case 2:{
			dir = "hello";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HelloContainer(i);
			}			
			break;
		}
		case 3:{
			dir = "searchlight";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new SearchlightContainer(i);
			}
			break;
		}
		case 4:{
			dir = "uconnect";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new UconnectContainer(i);
			}
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
	
	initEnergySupply2Star(theta0, dir); 
	
	stringstream ss;
	ss << dir<<"/star_energy2.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	w[0]->setTheta(theta0);
	w[0]->resetStartTime();
	for(int k = 0; k < TIMES && w[0]->alive == true; k++){
		bound = k;
		for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval && w[0]->alive == true; timeFrame ++){
			w[0]->endTime ++;
			if(w[0]->isOff(timeFrame) == false  && w[0]->alive == true){
				w[0]->energy -= pn;
				if(w[0]->energy <= 0 && flag == true){
					w[0]->alive = false;
					flag = false;
					break;
				}
				for(int j = 1; j < NODENUM; j++){
					if(w[j]->isOff(timeFrame) == false && w[j]->alive == true){
						onNeighbor++;
					}
				}
				if(onNeighbor == 1){
					for(int j = 1; j < NODENUM; j++){
						if(w[j]->isOff(timeFrame) == false && findNeighborTime[j] == 0 && w[j]->alive == true){
							findNeighborTime[j] = timeFrame;
							max = timeFrame;
						}
					}
				}
			}
			onNeighbor = 0;
		}
		
		w[0]->energySupply2();
		dutycycle[k] += w[0]->theta;
		percent[k] += w[0]->energy*1.0 / Pmax;
		for(int l = 1; l < NODENUM; l++){
			findNeighborTime[l] = 0;
		}
		if(max == 0){
			late[k] += timeInterval - w[0]->startTime;
		}
		else{
			late[k] += max - w[0]->startTime;
		}
		max = 0;
	}
	for(int k = 0; k < bound; k++){
		f<<dutycycle[k]<<","<<percent[k]<<","<<late[k]<<endl;
		
	}
	f<<"life Cycle = "<<w[0]->endTime - w[0]->startTime<<endl;
}



/*  Input: timeInterval: the number of time slots in a time period
 *		  alg: which algorithm
 *	Output: duty cycle, remaining energy percentage, average latency after each timeInterval
 *	Function: calculate the average number for HOW times of energy0 in star networks
 */
void multi0(int timeInterval, int alg){
	string dir;
	switch(alg){
		case 1:{
			dir = "hedis";
			break;
		}
		case 2:{
			dir = "hello";
			break;
		}
		case 3:{
			dir = "searchlight";
			break;
		}
		case 4:{
			dir = "uconnect";
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
	stringstream ss;
	ss <<dir<<"/star_ave_energy0.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	int HOW = 20;
	
	for(int k = 0; k < TIMES; k++){
		dutycycle[k] = 0;
		percent[k] = 0;
		late[k] = 0;
	}
	
	for(int k = 0; k < HOW; k++){
		energy_0(timeInterval, alg);
		
	}
	for(int k = 0; k < 20; k++){
		cout<<"k="<<k<<" "<<dutycycle[k]/HOW<<","<<percent[k]/HOW<<","<<late[k]/HOW<<endl;
		f<<dutycycle[k]/HOW<<","<<percent[k]/HOW<<","<<late[k]/HOW<<endl;
		
	}
	f<<"life Cycle = "<<w[0]->endTime - w[0]->startTime<<endl;
}


/*  Input: timeInterval: the number of time slots in a time period
 *		  alg: which algorithm
 *	Output: duty cycle, remaining energy percentage, average latency after each timeInterval
 *	Function: calculate the average number for HOW times of energy1 in star networks
 */
void multi1(int timeInterval, int alg){
	string dir;
	switch(alg){
		case 1:{
			dir = "hedis";
			break;
		}
		case 2:{
			dir = "hello";
			break;
		}
		case 3:{
			dir = "searchlight";
			break;
		}
		case 4:{
			dir = "uconnect";
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
	
	stringstream ss;
	ss <<dir<< "/star_ave_energy1.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	int HOW = 20;
	
	for(int k = 0; k < TIMES; k++){
		dutycycle[k] = 0;
		percent[k] = 0;
		late[k] = 0;
	}
	
	for(int k = 0; k < HOW; k++){
		energy_1(timeInterval, alg);
		
	}
	for(int k = 0; k < 20; k++){
		cout<<"k="<<k<<" "<<dutycycle[k]/HOW<<","<<percent[k]/HOW<<","<<late[k]/HOW<<endl;
		f<<dutycycle[k]/HOW<<","<<percent[k]/HOW<<","<<late[k]/HOW<<endl;
		
	}
	f<<"life Cycle = "<<w[0]->endTime - w[0]->startTime<<endl;
} 

/*  Input: timeInterval: the number of time slots in a time period
 *		  alg: which algorithm
 *	Output: duty cycle, remaining energy percentage, average latency after each timeInterval
 *	Function: calculate the average number for HOW times of energy2 in star networks
 */
void multi2(int timeInterval, int alg){
	string dir;
	switch(alg){
		case 1:{
			dir = "hedis";
			break;
		}
		case 2:{
			dir = "hello";
			break;
		}
		case 3:{
			dir = "searchlight";
			break;
		}
		case 4:{
			dir = "uconnect";
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return;
	}
	stringstream ss;
	ss <<dir<< "/star_ave_energy2.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	int HOW = 20;
	
	for(int k = 0; k < TIMES; k++){
		dutycycle[k] = 0;
		percent[k] = 0;
		late[k] = 0;
	}
	
	for(int k = 0; k < HOW; k++){
		energy_2(timeInterval, alg);
		
	}
	for(int k = 0; k < 20; k++){
		cout<<"k="<<k<<" "<<dutycycle[k]/HOW<<","<<percent[k]/HOW<<","<<late[k]/HOW<<endl;
		f<<dutycycle[k]/HOW<<","<<percent[k]/HOW<<","<<late[k]/HOW<<endl;
		
	}
	f<<"life Cycle = "<<w[0]->endTime - w[0]->startTime<<endl;
} 


/*  Input: vs: source node identifier
 *	Output: prev[i]: the previous node to go to node i
 *          dist[i]: the smallest distance between node vs and node i
 *	Function: calculate the smallest path and distance from source node vs 
 * to all the other nodes in throughput network with NODENUMPUT*NODENUMPUT nodes
 */
void Dijkstra(int vs, int prev[], int dist[])
{
    int i,j,k;
    int min;
    int tmp;
    int flag[NODENUMPUT];      
    
    for (i = 0; i < NODENUMPUT; i++)
    {
        flag[i] = 0;              
        prev[i] = 0;             
        dist[i] = (latePUT[vs][i] == INF) ? INF : 1;
    }

    flag[vs] = 1;
    dist[vs] = 0;
	k=-1;
    for (i = 1; i < NODENUMPUT; i++)
    {
        min = INF;
        for (j = 0; j < NODENUMPUT; j++)
        {
            if (flag[j]==0 && dist[j]<min)
            {
                min = dist[j];
                k = j;
            }
        }
        if(k==-1) break;
        flag[k] = 1;

        for (j = 0; j < NODENUMPUT; j++)
        {
            tmp = (latePUT[k][j]==INF ? INF : (min + 1));
            if (flag[j] == 0 && (tmp  < dist[j]) )
            {
                dist[j] = tmp;
                prev[j] = k;
            }
        }
    }
}


/*  Input: timeInterval: the number of time slots in a time period
 *         arrayLength: times to run and calculate the average
 *		   dir: directory to store output files
 *	Output: line i: the throughput at time slot i
 *	Function: calculate the throughput of bare protocol of 
 *  NODENUM*NODENUM aggregated network 
 */
void throughtput0(int timeInterval, int arrayLength, string dir)
{
	int onNeighbor[NODENUMPUT];
	stringstream ss;
	ss <<dir<< "/throughput0.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	initLatePUT();
	for(int i = 0; i < NODENUMPUT; i++){
		warray[i]->reset();
		warray[i]->setTheta(theta0);
		onNeighbor[i] = 0;
	}
	
	for(int k = 0; k < arrayLength; k++){
		cout<<"*********"<<k<<"***********"<<endl;
		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame++){
			for(int i = 0; i < NODENUMPUT; i++){
				warray[i]->endTime++;
				if(warray[i]->isOff(timeFrame) == false && warray[i]->alive == true){
					warray[i]->energy -= pn;
					if(warray[i]->energy <= 0){
						warray[i]->alive = false;
					}
					for(int j = 0; j < NODENUMPUT; j++){
						if(warray[j]->isOff(timeFrame) == false && neighborMapPUT[i][j] == 1  && warray[j]->alive == true){
							onNeighbor[j]++;
						}
					}
				}
			}

			for(int i = 0; i < NODENUMPUT; i++){
				if(onNeighbor[i] == 1 && warray[i]->isOff(timeFrame) ==false && warray[i]->alive == true){
					for(int j = 0; j < NODENUMPUT; j++){
						if(neighborMapPUT[i][j] == 1 && onNeighbor[j] == 1 && warray[j]->isOff(timeFrame) ==false && warray[j]->alive == true){
							setLatePUT(i, j, timeFrame);
						}
					}
				}
			}
			for(int i = 0; i < NODENUMPUT; i++){
				onNeighbor[i] = 0;
			}
			double result = 0.0;
			for(int i = 0; i < NODENUMPUT; i++){
				if(warray[i]->alive == true){
					int prev[NODENUMPUT] = {0};
	    			int dist[NODENUMPUT] = {0};
					Dijkstra(i, prev, dist);
					for(int j = 0; j < NODENUMPUT; j++){
						if(dist[j] < INF && dist[j] != 0 && warray[j]->alive == true){
							result=result+500.0/dist[j];
						}
					}
				}
			}
			f<<result<<endl;
		}
	}
}


/*  Input: timeInterval: the number of time slots in a time period
 *         arrayLength: times to run and calculate the average
 *		   dir: directory to store output files
 *	Output: line i: the throughput at time slot i
 *	Function: calculate the throughput of PDR of 
 *  NODENUM*NODENUM aggregated network 
 */
void throughtput1(int timeInterval, int arrayLength, string dir)
{
	int onNeighbor[NODENUMPUT];
	stringstream ss;
	ss <<dir<< "/throughput1.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);


	initLatePUT();
	for(int i = 0; i < NODENUMPUT; i++){
		warray[i]->reset();
		warray[i]->setTheta(theta0);
		onNeighbor[i] = 0;
	}

	
	for(int k = 0; k < arrayLength; k++){
		cout<<"*********"<<k<<"***********"<<endl;
		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame++){
			for(int i = 0; i < NODENUMPUT; i++){
				warray[i]->endTime++;
				if(warray[i]->isOff(timeFrame) == false&& warray[i]->alive == true){
					warray[i]->energy -= pn; 
					for(int j = 0; j < NODENUMPUT; j++){
						if(warray[j]->isOff(timeFrame) == false && neighborMapPUT[i][j] == 1&& warray[j]->alive == true){
							onNeighbor[j]++;
						}
					}
				}
			}

			for(int i = 0; i < NODENUMPUT; i++){
				if(onNeighbor[i] == 1 && warray[i]->isOff(timeFrame) == false&& warray[i]->alive == true){
					for(int j = 0; j < NODENUMPUT; j++){
						if(neighborMapPUT[i][j] == 1 && onNeighbor[j] == 1 && warray[j]->isOff(timeFrame) ==false && latePUT[i][j] == INF && warray[j]->alive == true){
							setLatePUT(i, j, timeFrame);
						}
					}
				}
			}
			for(int i = 0; i < NODENUMPUT; i++){
				onNeighbor[i] = 0;
			}
			double result = 0.0;
			for(int i = 0; i < NODENUMPUT; i++){
				if(warray[i]->alive == true){
					int prev[NODENUMPUT] = {0};
	    			int dist[NODENUMPUT] = {0};
					Dijkstra(i, prev, dist);
					for(int j = 0; j < NODENUMPUT; j++){
						if(dist[j] < INF && dist[j] != 0 && warray[j]->alive == true){
							result=result+500.0/dist[j];
						}
					}
				}
			}
			f<<result<<endl;
		}
		
		for(int i = 0; i < NODENUMPUT; i++){
			warray[i]->energySupply1(theta0);
		}
	}
}



/*  Input: timeInterval: the number of time slots in a time period
 *         arrayLength: times to run and calculate the average
 *		   dir: directory to store output files
 *	Output: line i: the throughput at time slot i
 *	Function: calculate the throughput of PWR of 
 *  NODENUM*NODENUM aggregated network 
 */
void throughtput2(int timeInterval, int arrayLength, string dir)
{
	int onNeighbor[NODENUMPUT];
	stringstream ss;
	ss <<dir<< "/throughput2.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);


	initLatePUT();
	for(int i = 0; i < NODENUMPUT; i++){
		warray[i]->reset();
		warray[i]->setTheta(theta0);
		onNeighbor[i] = 0;
	}
    initEnergySupply2PUT(theta0, dir);

	for(int k = 0; k < arrayLength; k++){
		cout<<"*********"<<k<<"***********"<<endl;
		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame++){
			for(int i = 0; i < NODENUMPUT; i++){
				warray[i]->endTime++;
				if(warray[i]->isOff(timeFrame) == false&& warray[i]->alive == true){
					warray[i]->energy -= pn;
					for(int j = 0; j < NODENUMPUT; j++){
						if(warray[j]->isOff(timeFrame) == false && neighborMapPUT[i][j] == 1&& warray[j]->alive == true){
							onNeighbor[j]++;
						}
					}
				}
			}

			for(int i = 0; i < NODENUMPUT; i++){
				if(onNeighbor[i] == 1 && warray[i]->isOff(timeFrame) == false&& warray[i]->alive == true){
					for(int j = 0; j < NODENUMPUT; j++){
						if(neighborMapPUT[i][j] == 1 && onNeighbor[j] == 1 && warray[j]->isOff(timeFrame) ==false && latePUT[i][j] == INF&& warray[j]->alive == true){
							setLatePUT(i, j, timeFrame);
						}
					}
				}
			}
			for(int i = 0; i < NODENUMPUT; i++){
				onNeighbor[i] = 0;
			}
			double result = 0.0;
			for(int i = 0; i < NODENUMPUT; i++){
				if(warray[i]->alive == true){
					int prev[NODENUMPUT] = {0};
	    			int dist[NODENUMPUT] = {0};
					Dijkstra(i, prev, dist);
					for(int j = 0; j < NODENUMPUT; j++){
						if(dist[j] < INF && dist[j] != 0 && warray[j]->alive == true){
							result=result+500.0/dist[j];
						}
					}
				}
			}
			f<<result<<endl;
		}
		
		for(int i = 0; i < NODENUMPUT; i++){
			warray[i]->energySupply2();
		}
	}
}

/*  Input: timeInterval: the number of time slots in a time period
 *         arrayLength: times to run and calculate the average
 *		   alg: which algorithm
 *	Output: 3 output files, bare protocol, PDR and PWR
 *	Function: output the throughput in 3 situations of
 *  NODENUM*NODENUM aggregated network 
 */
void throughput(int timeInterval, int arrayLength, int alg)
{
	string dir = initMapPUT(alg);
	if(dir == "") return;
	throughtput0(timeInterval, arrayLength, dir);
	throughtput1(timeInterval, arrayLength, dir);
	throughtput2(timeInterval, arrayLength, dir);
	
}


/* Input: prob: probability
 * Function: set warray's PROB1 to be prob
 */
void setwarrayProb1(double prob)
{
	for(int i = 0; i < 1000; i++){
		warray[i]->PROB1 = prob;
	}
}


/* Input: prob: probability
 * Function: set warray's PROB2 to be prob
 */
void setwarrayProb2(double prob)
{
	for(int i = 0; i < 1000; i++){
		warray[i]->PROB2 = prob;
	}
}



/* Input: timeInterval: the number of time slots in a time period
 *        arrayLength: times to run and calculate the average
 *		  dir: directory to store output files
 * Output: latency matrix
 * Function: calculate the latency matrix of bare protocol in 
 * 1000*1000 network 
 */
void distributionMod0(int timeInterval, int arrayLength, string dir)
{
	int onNeighbor[1000];
	stringstream ss;
	ss <<dir<< "/distribution_mod0_latency.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	initLatency();
    

    for(int i = 0; i < 1000; i++){
		onNeighbor[i] = 0;
		for(int j = 0; j < 1000; j++){
			aveLatency[i][j] = 0;
		}
	}
	for(int k = 0; k < arrayLength; k++){
		cout<<"******************"<<k<<"***************"<<endl;
		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame++){
			for(int i = 0; i < 1000; i++){
				if(warray[i]->isOff(timeFrame) == false){
					for(int j = 0; j < 1000; j++){
						if(warray[j]->isOff(timeFrame) == false && neighborMap[i][j] == 1){
							onNeighbor[j]++;
						}
					}
				}
			}

			for(int i = 0; i < 1000; i++){
				if(onNeighbor[i] == 1 && warray[i]->isOff(timeFrame) ==false){
					for(int j = 0; j < 1000; j++){
						if(neighborMap[i][j] == 1 && onNeighbor[j] == 1 && warray[j]->isOff(timeFrame) ==false){
							setLatency(i, j, timeFrame);
							setLatency(j, i, timeFrame);
						}
					}
				}
			}
			for(int i = 0; i < 1000; i++){
				onNeighbor[i] = 0;
			}
		}

		for(int i = 0; i < 1000; i++){
			for(int j = 0; j < 1000; j++){
				if(neighborMap[i][j] == 1 && latencyMatrix[i][j] != MAX_INT){
					aveLatency[i][j] += latencyMatrix[i][j];
				}
			}
		}
		initLatency();
	}
	for(int i = 0; i < 1000; i++){
		for(int j = 0; j < 1000; j++){
			aveLatency[i][j] = aveLatency[i][j]*1.0/arrayLength;
			f<<aveLatency[i][j]<<",";
		}
		f<<endl;
	}
}


/* Input: timeInterval: the number of time slots in a time period
 *        arrayLength: times to run and calculate the average
 *		  dir: directory to store output files
 * Output: latency matrix
 * Function: calculate the latency matrix of PPR in 
 * 1000*1000 network 
 */
void distributionMod1(int timeInterval, int arrayLength, string dir)
{
	int onNeighbor[1000];
	stringstream ss;
	ss <<dir<< "/distribution_mod1_latency.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);


	initLatency();

    for(int i = 0; i < 1000; i++){
		onNeighbor[i] = 0;
		for(int j = 0; j < 1000; j++){
			aveLatency[i][j] = 0;
		}
	}
	for(int k = 0; k < arrayLength; k++){
		cout<<"******************"<<k<<"***************"<<endl;
		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame++){

			for(int i = 0; i < 1000; i++){
				if(warray[i]->isOff_collision1(timeFrame) == false){
					for(int j = 0; j < 1000; j++){
						if(warray[j]->isOff_collision1(timeFrame) == false && neighborMap[i][j] == 1){
							onNeighbor[j]++;
						}
					}
				}
			}

			for(int i = 0; i < 1000; i++){
				if(onNeighbor[i] == 1 && warray[i]->isOff_collision1(timeFrame) ==false){
					for(int j = 0; j < 1000; j++){
						if(neighborMap[i][j] == 1 && onNeighbor[j] == 1 && warray[j]->isOff_collision1(timeFrame) ==false){
							setLatency(i, j, timeFrame);
							setLatency(j, i, timeFrame);
						}
					}
				}
			}
			for(int i = 0; i < 1000; i++){
				onNeighbor[i] = 0;
			}
		}

		for(int i = 0; i < 1000; i++){
			for(int j = 0; j < 1000; j++){
				if(neighborMap[i][j] == 1 && latencyMatrix[i][j] != MAX_INT){
					aveLatency[i][j] += latencyMatrix[i][j];
				}
			}
		}
		initLatency();
	}
	for(int i = 0; i < 1000; i++){
		for(int j = 0; j < 1000; j++){
			aveLatency[i][j] = aveLatency[i][j]*1.0/arrayLength;
			f<<aveLatency[i][j]<<",";
		}
		f<<endl;
	}
}



/* Input: timeInterval: the number of time slots in a time period
 *        arrayLength: times to run and calculate the average
 *		  dir: directory to store output files
 * Output: latency matrix
 * Function: calculate the latency matrix of DPR in 
 * 1000*1000 network 
 */
void distributionMod2(int timeInterval, int arrayLength, string dir)
{
	int onNeighbor[1000];
	stringstream ss;
	ss <<dir<< "/distribution_mod2_latency.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);

	initLatency();

    for(int i = 0; i < 1000; i++){
		onNeighbor[i] = 0;
		for(int j = 0; j < 1000; j++){
			aveLatency[i][j] = 0;
		}
	}
	for(int k = 0; k < arrayLength; k++){
		cout<<"******************"<<k<<"***************"<<endl;
		for(int timeFrame = 1; timeFrame <= timeInterval; timeFrame++){
			for(int i = 0; i < 1000; i++){
				if(warray[i]->isOff_collision2(timeFrame, timeInterval) == false){
					for(int j = 0; j < 1000; j++){
						if(warray[j]->isOff_collision2(timeFrame, timeInterval) == false && neighborMap[i][j] == 1){
							onNeighbor[j]++;
						}
					}
				}
			}

			for(int i = 0; i < 1000; i++){
				if(onNeighbor[i] == 1 && warray[i]->isOff_collision2(timeFrame, timeInterval) ==false){
					for(int j = 0; j < 1000; j++){
						if(neighborMap[i][j] == 1 && onNeighbor[j] == 1 && warray[j]->isOff_collision2(timeFrame, timeInterval) ==false){
							setLatency(i, j, timeFrame);
							setLatency(j, i, timeFrame);
						}
					}
				}
			}
			for(int i = 0; i < 1000; i++){
				onNeighbor[i] = 0;
			}
		}

		for(int i = 0; i < 1000; i++){
			for(int j = 0; j < 1000; j++){
				if(neighborMap[i][j] == 1 && latencyMatrix[i][j] != MAX_INT){
					aveLatency[i][j] += latencyMatrix[i][j];
				}
			}
		}
		initLatency();
	}
	for(int i = 0; i < 1000; i++){
		for(int j = 0; j < 1000; j++){
			aveLatency[i][j] = aveLatency[i][j]*1.0/arrayLength;
			f<<aveLatency[i][j]<<",";
		}
		f<<endl;
	}
}


/* Input: timeInterval: the number of time slots in a time period
 *        arrayLength: times to run and calculate the average
 *		  alg: which algorithm
 * Output: 3 latency matrix
 * Function: calculate the latency matrix of bare protocol, PPR, DPR in 
 * 1000*1000 network 
 */
void distribution(int timeInterval, int arrayLength, int alg)
{
	string dir = initMapDistr(alg);
	if(dir == "") return;
	setwarrayProb1(0.4);
	setwarrayProb2(0.2);
	distributionMod0(timeInterval, arrayLength, dir);
	distributionMod1(timeInterval, arrayLength, dir);
	distributionMod2(timeInterval, arrayLength, dir);
}



/* Input: alg: which algorithm
 * Output: dir: directory to store output files
 * Function: get responding dir and wificontainer from alg in star networks 
 */
string initw(int alg)
{
	string dir="";
	switch(alg){
		case 1:{
			dir = "hedis";
			for(int i = 0; i < NODENUMBOUND; i++){
				w[i] = new HedisContainer(i);
			}
			return dir;
		}
		case 2:{
			dir = "hello";
			for(int i = 0; i < NODENUMBOUND; i++){
				w[i] = new HelloContainer(i);
			}			
			return dir;
		}
		case 3:{
			dir = "searchlight";
			for(int i = 0; i < NODENUMBOUND; i++){
				w[i] = new SearchlightContainer(i);
			}
			return dir;
		}
		case 4:{
			dir = "uconnect";
			for(int i = 0; i < NODENUMBOUND; i++){
				w[i] = new UconnectContainer(i);
			}
			return dir;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			return dir;
	}
	
}


/* Input: prob: probability
 * Function: set w's PROB1 to be prob
 */
void setwProb1(double prob)
{
	for(int i = 0; i < NODENUMBOUND; i++){
		w[i]->PROB1 = prob;
	}
}



/* Input: prob: probability
 * Function: set w's PROB2 to be prob
 */
void setwProb2(double prob)
{
	for(int i = 0; i < NODENUMBOUND; i++){
		w[i]->PROB2 = prob;
	}
}



/* Input: timeInterval: the number of time slots in a time period
 *        dir: directory to store output files
 * Output: the percentage of node that could be found as node number changes
 * Function: calculate the discovery rate of bare protocol in star networks 
 */
void nodenumMod0(int timeInterval, string dir)
{
	int findNeighborTime[NODENUMBOUND];
	int onNeighbor = 0;
	int count = 0; 
	
	for(int i = 0; i < NODENUMBOUND; i++){
		findNeighborTime[i] = 0;
	}
	
	stringstream ss;
	ss <<dir<< "/nodenum_mod0.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	w[0]->setTheta(theta0);
	w[0]->resetStartTime();
	
	for(int nodenum = 2; nodenum < NODENUMBOUND; nodenum ++){
		for(int k = 0; k < TIMES; k++){
			for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval; timeFrame ++){
				if(w[0]->isOff(timeFrame) == false){
					for(int j = 1; j < nodenum; j++){
						if(w[j]->isOff(timeFrame) == false && w[j]->alive == true){
							onNeighbor++;
						}
					}
					if(onNeighbor == 1){
						for(int j = 1; j < nodenum; j++){
							if(w[j]->isOff(timeFrame) == false && findNeighborTime[j] == 0 && w[j]->alive == true){
								findNeighborTime[j] = timeFrame;
							}
						}
					}
				}
				onNeighbor = 0;
			}
			for(int j = 1; j < nodenum; j++){
				if(findNeighborTime[j] != 0){
					count++;
				}
			}
			for(int l = 1; l < nodenum; l++){
				findNeighborTime[l] = 0;
			}
			
		}
		f<<nodenum-1<<","<<count*1.0/TIMES<<endl;
		count = 0;
	}
}



/* Input: timeInterval: the number of time slots in a time period
 *        dir: directory to store output files
 * Output: the percentage of node that could be found as node number changes
 * Function: calculate the discovery rate of PPR in star networks 
 */
void nodenumMod1(int timeInterval, string dir)
{
	int findNeighborTime[NODENUMBOUND];
	int onNeighbor = 0;
	int count = 0; 
	
	for(int i = 0; i < NODENUMBOUND; i++){
		findNeighborTime[i] = 0;
	}
	
	stringstream ss;
	ss <<dir<< "/nodenum_mod1.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	w[0]->setTheta(theta0);
	w[0]->resetStartTime();
	
	for(int nodenum = 2; nodenum < NODENUMBOUND; nodenum ++){
		for(int k = 0; k < TIMES; k++){
			for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval; timeFrame ++){
				if(w[0]->isOff_collision1(timeFrame) == false){
					for(int j = 1; j < nodenum; j++){
						if(w[j]->isOff_collision1(timeFrame) == false && w[j]->alive == true){
							onNeighbor++;
						}
					}
					if(onNeighbor == 1){
						for(int j = 1; j < nodenum; j++){
							if(w[j]->isOff_collision1(timeFrame) == false && findNeighborTime[j] == 0 && w[j]->alive == true){
								findNeighborTime[j] = timeFrame;
							}
						}
					}
				}
				onNeighbor = 0;
			}
			for(int j = 1; j < nodenum; j++){
				if(findNeighborTime[j] != 0){
					count++;
				}
			}
			for(int l = 1; l < nodenum; l++){
				findNeighborTime[l] = 0;
			}
			
		}
		f<<nodenum-1<<","<<count*1.0/TIMES<<endl;
		count = 0;
	}
}



/* Input: timeInterval: the number of time slots in a time period
 *        dir: directory to store output files
 * Output: the percentage of node that could be found as node number changes
 * Function: calculate the discovery rate of DPR in star networks 
 */
void nodenumMod2(int timeInterval, string dir)
{
	int findNeighborTime[NODENUMBOUND];
	int onNeighbor = 0;
	int count = 0;
	bool temp[NODENUMBOUND];
	
	for(int i = 0; i < NODENUMBOUND; i++){
		findNeighborTime[i] = 0;
	}
	
	stringstream ss;
	ss <<dir<< "/nodenum_mod2.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	
	w[0]->setTheta(theta0);
	w[0]->resetStartTime();
	for(int nodenum = 2; nodenum < NODENUMBOUND; nodenum ++){
		for(int k = 0; k < TIMES; k++){
			for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval; timeFrame ++){
				for(int x = 0; x < nodenum; x++){
					temp[x] = w[x]->isOff_collision2(timeFrame, timeInterval);
				}
				if(temp[0] == false  && w[0]->alive == true){
					for(int j = 1; j < nodenum; j++){
						if(temp[j] == false && w[j]->alive == true){
							onNeighbor++;
						}
					}
					if(onNeighbor == 1){
						for(int j = 1; j < nodenum; j++){
							if(temp[j] == false && findNeighborTime[j] == 0 && w[j]->alive == true){
								findNeighborTime[j] = timeFrame;
							}
						}
					}
				}
				onNeighbor = 0;
			}
			for(int j = 1; j < nodenum; j++){
				if(findNeighborTime[j] != 0){
					count++;
				}
			}
			for(int l = 1; l < nodenum; l++){
				findNeighborTime[l] = 0;
			}
		}
		cout<<nodenum-1<<","<<count*1.0/TIMES<<endl;
		
		f<<nodenum-1<<","<<count*1.0/TIMES<<endl;
		count = 0;
	}

}



/* Input: timeInterval: the number of time slots in a time period
 *        alg: which algorithm
 * Output: discovery rate of 3 situations
 * Function: calculate the discovery rate of bare protocol, PPR and DPR in star networks 
 */
void nodenum(int timeInterval, int alg)
{
	string dir = initw(alg);
	if(dir == "") return;
	setwProb1(0.4);
	setwProb2(0.2);
	nodenumMod0(timeInterval, dir);
	nodenumMod1(timeInterval, dir);
	nodenumMod2(timeInterval, dir);
}



/* Input: the: duty cycle
 * Function: set duty cycle of w to be the
 */
void setAllTheta(double the){
	for(int i = 0; i < NODENUM; i++){
		w[i]->setTheta(the);
	}
}



/* Input: timeInterval: the number of time slots in a time period
 *        thetaaa: the duty cycle of all nodes
 *        dir: directory to store output files
 * Output: average latency of different PROB1
 * Function: calculate the average latency with PROB1's changes of PPR in star networks 
 */
void sensitivityP1(int timeInterval, double thetaaa, string dir)
{
	int findNeighborTime[NODENUM];
	int starlatency = MAX_INT;;
	int onNeighbor = 0;
  	int sum = 0;
  	
	for(int i = 0; i < NODENUM; i++){
		findNeighborTime[i] = 0;
	}
	stringstream ss;
	ss <<dir<< "/sensitivity_p1.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	f<<thetaaa<<endl;
	setAllTheta(thetaaa);
	for(double prob = 0.1; prob <= 0.9; prob+=0.1){
		setwProb1(prob);
		for(int s = 0; s < TIMES; s++){
			for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval; timeFrame ++){
				if(w[0]->isOff_collision1(timeFrame) == false){
					for(int j = 1; j < NODENUM; j++){
						if(w[j]->isOff_collision1(timeFrame) == false){
							onNeighbor++;
						}
					}
					if(onNeighbor == 1){
						for(int j = 1; j < NODENUM; j++){
							if(w[j]->isOff_collision1(timeFrame) == false && findNeighborTime[j] == 0){
								findNeighborTime[j] = timeFrame;
							}
						}
					}
				}
				onNeighbor = 0;
			}

		    for(int l = 1; l < NODENUM; l++){
		          if(findNeighborTime[l] == 0){
		            sum+=timeInterval;
		          }
		          else{
		            sum+=findNeighborTime[l];
		            findNeighborTime[l] = 0;
		          }
		
		    }

	    }
	    f<<prob<<","<<sum*1.0/NODENUM/TIMES<<endl;
	    sum=0;
	}
}


/* Input: timeInterval: the number of time slots in a time period
 *        thetaaa: the duty cycle of all nodes
 *        dir: directory to store output files
 * Output: average latency of different PROB1
 * Function: calculate the average latency with PROB2's changes of DPR in star networks 
 */
void sensitivityP2(int timeInterval, double thetaaa, string dir)
{
	int findNeighborTime[NODENUM];
	int onNeighbor = 0;
	int sum=0;
	bool temp[NODENUM];
	for(int i = 0; i < NODENUM; i++){
		findNeighborTime[i] = 0;
	}
	
	stringstream ss;
	ss <<dir<< "/sensitivity_p2.txt";
	fstream f (ss.str().c_str(), std::ios::out | std::ios::app);
	f<<thetaaa<<endl;
	setAllTheta(thetaaa);
	
	for(double prob = 0.1; prob <= 0.9; prob+=0.1) {
		setwProb2(prob);
		for(int s = 0; s < TIMES; s++){
			for(int timeFrame = w[0]->startTime; timeFrame <= timeInterval; timeFrame ++){
				for(int x = 0; x < NODENUM; x++){
					temp[x] = w[x]->isOff_collision2(timeFrame, timeInterval);
				}
				if(temp[0] == false){
					for(int j = 1; j < NODENUM; j++){
						if(temp[j] == false){
							onNeighbor++;
						}
					}
					if(onNeighbor == 1){
						for(int j = 1; j < NODENUM; j++){
							if(temp[j] == false && findNeighborTime[j] == 0){
								findNeighborTime[j] = timeFrame;
							}
						}
					}
				}
				onNeighbor = 0;
			}
			
			for(int l = 1; l < NODENUM; l++){
				if(findNeighborTime[l]==0){
					sum+=timeInterval;
				}
				else{
					sum+=findNeighborTime[l];
					findNeighborTime[l] = 0;
				}
			}
			for(int l = 1; l < NODENUM; l++){
				findNeighborTime[l] = 0;
			}
		}
		f<<prob<<","<<sum/TIMES/(NODENUM-1);
		sum=0;
	}
}



/* Input: timeInterval: the number of time slots in a time period
 *        thetaaa: the duty cycle of all nodes
 *        alg: which algorithm
 * Output: average latency of an algrithm with different PROB1, PROB2
 * Function: sensitivity analysis about PROB1 of PPR and PROB2 in DPR
 */
void sensitivity(int timeInterval, double thetaaa, int alg)
{
	string dir="";
	switch(alg){
		case 1:{
			dir = "hedis";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HedisContainer(i);
			}
			break;
		}
		case 2:{
			dir = "hello";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new HelloContainer(i);
			}			
			break;
		}
		case 3:{
			dir = "searchlight";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new SearchlightContainer(i);
			}
			break;
		}
		case 4:{
			dir = "uconnect";
			for(int i = 0; i < NODENUM; i++){
				w[i] = new UconnectContainer(i);
			}
			break;
		}
		default:
			cout<<"Error: no such algorithms!"<<endl;
			break;
	}
	if(dir == "") return;
	sensitivityP1(timeInterval, thetaaa, dir);
	sensitivityP2(timeInterval, thetaaa, dir);
}



int main(int argc, char** argv) {
	int alg;
	int net;
	double thetaaa;
	
	srand(time(NULL));
	findAllPrimes(10000);

	
	while(true){
		cout<<"Welcome to Spear!\nHere is the simulation part. "
		"Spear can evaluate the latency of two nodes:\n"
		"1 two Node(symmetric)\n2 two Node(asymmetric)\n"
		"can save power with PWR and PDR:\n"
		"3 energy and life time in star network\n4 throughput in "<<NODENUMPUT<<"*"<<NODENUMPUT<<" network\n"
		"can reduce collisions with PPR and DPR:\n"
		"5 discovery rate in 1000*1000 network\n"
		"6 discovery rate with different number of neighbors in star network\n"
		"7 sensitivity analysis in star network\n";
		cout<<"Please input what situation you want to run:\n";
		cin>>net;
		switch(net){
			case 1:{
				cout<<"Please input which algorithm you want to run:\n1 Hedis\n2 Hello\n3 Searchlight\n4 U-Connect\n5 Quorum\n6 all\n7 Exit"<<endl;
				cin>>alg;
				if(alg == 6){
					for(int i = 1; i < alg; i++){
						twoNodeSym(100000, 100, i);
					}
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 1 || alg == 2 || alg == 3 || alg == 4 || alg == 5){
					twoNodeSym(100000, 100, alg);
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 7){
					
				}
				cout<<"********return********"<<endl<<endl;
				break;
			}
			case 2:{
				cout<<"Please input which algorithm you want to run:\n1 Hedis\n2 Hello\n3 Searchlight\n4 U-Connect\n5 all\n6 Exit"<<endl;
				cin>>alg;
				if(alg == 5){
					for(int i = 1; i < alg; i++){
						twoNodeAsym(100000, 100, i);
					}
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 1 || alg == 2 || alg == 3 || alg == 4){
					twoNodeAsym(100000, 100, alg);
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 6){
				}
				cout<<"********return********"<<endl<<endl;
				break;
			}
			case 3:{
				cout<<"Please input which algorithm you want to run:\n1 Hedis\n2 Hello\n3 Searchlight\n4 U-Connect\n5 all\n6 Exit"<<endl;
				cin>>alg;
				if(alg == 5){
					for(int i = 1; i < alg; i++){
						multi0(100000, i);
						multi1(100000, i);
						multi2(100000, i);
					}	
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 1 || alg == 2 || alg == 3 || alg == 4){
					multi0(100000, alg);
					multi1(100000, alg);
					multi2(100000, alg);
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 6){
				}
				cout<<"********return********"<<endl<<endl;
				break;
			}
			case 4:{
				cout<<"Please input which algorithm you want to run:\n1 Hedis\n2 Hello\n3 Searchlight\n4 U-Connect\n5 all\n6 Exit"<<endl;
				cin>>alg;
				if(alg == 5){
					for(int i = 1; i < alg; i++){
						throughput(10000, 200, i);
					}
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 1 || alg == 2 || alg == 3 || alg == 4){
					throughput(10000, 200, alg);
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 6){
					
				}
				cout<<"********return********"<<endl<<endl;
				break;
			}
			case 5:{
				cout<<"Please input which algorithm you want to run:\n1 Hedis\n2 Hello\n3 Searchlight\n4 U-Connect\n5 all\n6 Exit"<<endl;
				cin>>alg;
				if(alg == 5){
					for(int i = 1; i < alg; i++){
						distribution(10000, 20, i);
					}
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 1 || alg == 2 || alg == 3 || alg == 4){
					distribution(10000, 20, alg);
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 6){
				}
				cout<<"********return********"<<endl<<endl;
				break;
			}
			case 6:{
				cout<<"Please input which algorithm you want to run:\n1 Hedis\n2 Hello\n3 Searchlight\n4 U-Connect\n5 all\n6 Exit"<<endl;
				cin>>alg;
				if(alg == 5){
					for(int i = 1; i < alg; i++){
						nodenum(10000, i);
					}
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 1 || alg == 2 || alg == 3 || alg == 4){
					nodenum(10000, alg);
					cout<<"Succeed!"<<endl;
				}
				else if(alg == 6){
				}
				cout<<"********return********"<<endl<<endl;
				break;
			}
			case 7:{
				cout<<"Please input which algorithm you want to run:\n1 Hedis\n2 Hello\n3 Searchlight\n4 U-Connect\n5 all\n6 Exit"<<endl;
				cin>>alg;
				cout<<"Please input your theta value:"<<endl;
				cin>>thetaaa;
				if(thetaaa > 0 && thetaaa < 1){
					if(alg == 5){
						for(int i = 1; i < alg; i++){
							sensitivity(100000, thetaaa, i);
						}
						cout<<"Succeed!"<<endl;
					}
					else if(alg == 1 || alg == 2 || alg == 3 || alg == 4){
						sensitivity(100000, thetaaa, alg);
						cout<<"Succeed!"<<endl;
					}
					else if(alg == 6){
					}
				}
				else{
					cout<<"Error: theta value!"<<endl;
				}
				cout<<"********return********"<<endl<<endl;
				break;
			}
			default:{
				cout<<"no such choice in the menu. Please try again.\n"; 
				cout<<"********return********"<<endl<<endl;
				break;
			}
		}
	
	}
	return 0;
}
