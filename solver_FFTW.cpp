/********************************
This is the CPP file of FFTW solver
**********************************/

#include "solver_FFTW.h"

Solver_FFTW::Solver_FFTW(){
	Input* initInput = new Input();
	numOfXGrid = initInput->getXGridNum();
	numOfYGrid = initInput->getYGridNum();

	
}

Solver_FFTW::~Solver_FFTW(){
}

void Solver_FFTW::firstDerivative(){

}

void Solver_FFTW::secondDerivative(){
}


// A function to calculate the energy
double Solver_FFTW::calculateE(){
	double E = 0;
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			E += 0.5*(ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j]);
		}
	}
	return E;

}

void Solver_FFTW::burgersSolver_FFTW(){
}