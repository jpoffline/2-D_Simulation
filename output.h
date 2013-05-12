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

/*===============================
Explanations to the flags(output type):

_VELOCITY: generating output for velocities.
_DERIVATIVEv: generating output for v_x and v_y.
_DERIVATIVEw: generating output for w_x and w_y.
_DDERIVATIVEv: generating output for v_x_x and v_y_y.
_DDERIVATIVEw: generating output for w_x_x and w_y_y.

===============================*/

enum output_Type{
	_VELOCITY,
	_DERIVATIVEv,
	_DERIVATIVEw,
	_DDERIVATIVEv,
	_DDERIVATIVEw
};

class Output{
public:
	Output(int numOfXGrid, int numOfYGrid, double **xVelocity, double **yVelocity, long int time,double energy,output_Type fileType);
	~Output();
};