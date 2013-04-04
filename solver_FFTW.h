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
	double **v;
	double **w;
	double **temp_Velocity;// temp is used in executing the plans
	double **firstD_u;
	double **secondD_u;

	//arrays to store first order derivative
	double **v_x;
	double **v_y;
	double **w_x;
	double **w_y;

	//arrays to store second order derivative
	double **v_x_x;
	double **v_y_y;
	double **w_x_x;
	double **w_y_y;

	//These are the arrays in fourier space
	fftw_complex **V;
	fftw_complex **W;
	fftw_complex **temp_U;
	fftw_complex **firstD_U;
	fftw_complex **secondD_U;

	//FFT transform
	fftw_plan plan_r2c;
	fftw_plan plan_c2r;
	fftw_plan plan_firstD;
	fftw_plan plan_secondD;

	double ***Adams_v;
	double ***Adams_w;
	ofstream energy;
};