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

class Solver_FFTW{
public:
	Solver_FFTW();
	~Solver_FFTW();
	void burgersSolver_FFTW();
private:
}£»