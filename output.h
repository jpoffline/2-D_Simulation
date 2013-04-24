/*********************************
This is the onput of the simulation.
The outputs files are named by the time and 
contain all the velocites of all girds.

Author: Francis Chen
Date: 04.03.2013
****************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "const.h"

using namespace std;

class Output{
public:
	Output(int numOfXGrid, int numOfYGrid, double **xVelocity, double **yVelocity, long int time);
	~Output();
};