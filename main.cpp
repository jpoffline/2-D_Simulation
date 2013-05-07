/**********************************
This is the main function of the simulation

Author: Francis Chen
Date: 04.03.2013
*********************************/

#include "solver_FFTW.h"

void main(){
	
	Solver_FFTW* test = new Solver_FFTW();
	test->burgersSolver_FFTW();
	return;
}