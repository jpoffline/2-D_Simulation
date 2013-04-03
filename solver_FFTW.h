/****************************************
This is the FFTW Solver for Burgers equation. FFTW is a package 
developed by MIT to solve problems by FFT.

Author: Francis Chen
Date: 04.03.2013
********************************************/

#include <fftw3.h>
#include "const.h"
#include "input.h"
#include "output.h"

using namespace std;

class Solver_FFTW{
public:
	Solver_FFTW();
	~Solver_FFTW();
	void burgersSolver_FFTW();
	void firstDerivative();
	void secondDerivative();
	double calculateE();

private:
	int numOfXGrid;
	int numOfYGrid;

	//These are the arrays in real space
	double **ux;
	double **uy;
	double **ux_x;
	double **ux_y;
	double **uy_x;
	double **uy_y;
	double **ux_x_x;
	double **ux_y_y;
	double **uy_x_x;
	double **ux_y_y;

	//These are the arrays in fourier space
	fftw_complex **U;
	fftw_complex **temp_U;
	fftw_complex **U_1;
	fftw_complex **U_2;

	//FFT transform
	fftw_plan plan_r2c;
	fftw_plan plan_c2r;
	fftw_plan plan_firstD;
	fftw_plan plan_secondD;

	double ***Adams;
	ofstream energy;
};