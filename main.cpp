/**********************************
This is the main function of the simulation
*********************************/

#include "solver_FFTW.h"

void main(){
	clock_t start = clock();
	Solver_FFTW* test = new Solver_FFTW();
	test->burgersSolver_FFTW();
	return;
}