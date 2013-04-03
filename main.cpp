/**********************************
This is the main function of the simulation
*********************************/

#include "solver_FFTW.h"

void main(){
	Solver_FFTW* test = new Solver_FFTW();
	test->burgersSolver_FFTW();
	return;
}