/*********************************
This is the input of the simulation.
In input, all the initial data will be read

Author: Francis Chen
Date: 04.03.2013
****************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include "const.h"

using namespace std;

class Input{
public:
	Input();
	~Input();
	int getXGridNum();
	int getYGridNum();
// i is the x coordinate and j is the y coordinate/
	double getVelocity(int i,int j); 

private:
	int numOfXGrid; //total number of grids in x axis
	int numOfYGrid; //total number of grids in y axis
	double **initVelocity; //array to store the input velocities

};