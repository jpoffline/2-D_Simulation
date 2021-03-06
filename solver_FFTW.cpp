/********************************
This is the CPP file of FFTW solver

Author: Francis Chen
Date: 04.03.2013
**********************************/

#include "solver_FFTW.h"


//constructor of the class
Solver_FFTW::Solver_FFTW(){
	/*===============================================*/
	Input* initInput = new Input();
	numOfXGrid = initInput->getXGridNum();
	numOfYGrid = initInput->getYGridNum();
	cout << "Input Finished" << endl;
	cout << "Grids along x axis: " << numOfXGrid << endl;
	cout << "Grids along y axis: " << numOfYGrid << endl;

	/*======================================================
	Initializing arrays in real space
	======================================================*/
	v = new double*[numOfXGrid];
	w = new double*[numOfXGrid];
	initE = 0;
	temp_Velocity = new double[numOfXGrid*numOfYGrid];
	firstD_u = new double[numOfXGrid*numOfYGrid];
	secondD_u = new double[numOfXGrid*numOfYGrid];
	for(int i = 0; i < numOfXGrid; i++){
		v[i] = new double[numOfYGrid];
		w[i] = new double[numOfYGrid];
		for(int j = 0; j < numOfYGrid; j++){
			v[i][j] = initInput->getXVelocity(i,j);
			w[i][j] = initInput->getYVelocity(i,j);
			initE += v[i][j]*v[i][j] + w[i][j]*w[i][j]; //calculating the initial energy
		}
	}

	for(int i = 0; i < numOfXGrid*numOfYGrid; i++){
		temp_Velocity[i] = 0;
		firstD_u[i] = 0;
		secondD_u[i] = 0;
	}
	
	/*=====================================================
	Initializing first order derivatives
	=====================================================*/
	v_x = new double*[numOfXGrid];
	v_y = new double*[numOfXGrid];
	w_x = new double*[numOfXGrid];
	w_y = new double*[numOfXGrid];
	for(int i = 0; i < numOfXGrid; i++){
		v_x[i] = new double[numOfYGrid];
		v_y[i] = new double[numOfYGrid];
		w_x[i] = new double[numOfYGrid];
		w_y[i] = new double[numOfYGrid];
		for(int j = 0; j < numOfYGrid; j++){
			v_x[i][j] = 0;
			v_y[i][j] = 0;
			w_x[i][j] = 0;
			w_y[i][j] = 0;
		}
	}
	
	/*=====================================================
	Initializing second order derivatives
	=====================================================*/
	v_x_x = new double*[numOfXGrid];
	v_y_y = new double*[numOfXGrid];
	w_x_x = new double*[numOfXGrid];
	w_y_y = new double*[numOfXGrid];
	for(int i = 0; i < numOfXGrid; i++){
		v_x_x[i] = new double[numOfYGrid];
		v_y_y[i] = new double[numOfYGrid];
		w_x_x[i] = new double[numOfYGrid];
		w_y_y[i] = new double[numOfYGrid];
		for(int j = 0; j < numOfYGrid; j++){
			v_x_x[i][j] = 0;
			v_y_y[i][j] = 0;
			w_x_x[i][j] = 0;
			w_y_y[i][j] = 0;
		}
	}

	/*========================================================
	Initializing the forces
	========================================================*/
	externalFx = new double*[numOfXGrid];
	externalFy = new double*[numOfXGrid];
	for(int i = 0;i < numOfXGrid; i++){
		externalFx[i] = new double[numOfYGrid];
		externalFy[i] = new double[numOfYGrid];
		for(int j = 0;j < numOfYGrid; j++){
			externalFx[i][j] = 0;
			externalFy[i][j] = 0;
		}
	}


	/*========================================================
	initializing multiple threads
	==========================================================*/
	if(fftw_init_threads()){
		fftw_plan_with_nthreads(THREADS);
		cout << "Using "<< THREADS << " threads" << endl << endl;
	}
	else {
		cout << "Using multiple threads failed" << endl;
		exit(0);
	}

	/*=====================================================
	Initializing arrays in fourier space
	=====================================================*/
	V = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*)*numOfXGrid);
	W = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*)*numOfXGrid);
	temp_U = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*numOfXGrid*(numOfYGrid/2+1));
	firstD_U = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*numOfXGrid*(numOfYGrid/2+1));
	secondD_U = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*numOfXGrid*(numOfYGrid/2+1));
	for(int i = 0; i < numOfXGrid; i++){
		V[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(numOfYGrid/2+1));
		W[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(numOfYGrid/2+1));
		for(int j = 0; j < numOfYGrid/2 + 1; j++){
			V[i][j][0] = 0;
			V[i][j][1] = 0;
			W[i][j][0] = 0;
			W[i][j][1] = 0;
		}
	}

	for(int i = 0; i < numOfXGrid*(numOfYGrid/2+1); i++){
		temp_U[i][0] = 0;
		temp_U[i][1] = 0;
		firstD_U[i][0] = 0;
		firstD_U[i][1] = 0;
		secondD_U[i][0] = 0;
		secondD_U[i][1] = 0;
	}

	temp = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*)*numOfXGrid);
	for(int i = 0; i < numOfXGrid; i++){
		temp[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(numOfYGrid/2+1));
	}
	/**=====================================================
	Initializing Plans
	======================================================*/
	plan_r2c = fftw_plan_dft_r2c_2d(numOfXGrid,numOfYGrid,temp_Velocity,temp_U,FFTW_ESTIMATE);
	plan_c2r = fftw_plan_dft_c2r_2d(numOfXGrid,numOfYGrid,temp_U,temp_Velocity,FFTW_ESTIMATE);
	plan_firstD = fftw_plan_dft_c2r_2d(numOfXGrid,numOfYGrid,firstD_U,firstD_u,FFTW_ESTIMATE);
	plan_secondD = fftw_plan_dft_c2r_2d(numOfXGrid,numOfYGrid,secondD_U,secondD_u,FFTW_ESTIMATE);

	/*======================================================
	initializing output energy file
	======================================================*/

	stringstream fileNameOfE;
	fileNameOfE << OUTPUT_PATH << "E" << ".txt";
	energy.open(fileNameOfE.str().c_str());

	/*======================================================
	generating the readMe.txt file
	======================================================*/

	stringstream fileNameOfReadMe;
	fileNameOfReadMe << OUTPUT_PATH << "readMe.txt";
	readMe.open(fileNameOfReadMe.str().c_str());
	readMe << "Nx = " << numOfXGrid << endl;
	readMe << "Ny = " << numOfYGrid << endl;
	readMe << "dt = " << TIME_STEP << endl;
	readMe << "Initial Energy = " << initE << endl;
	readMe << "Viscosity = "<< VISCOSITY << endl;
	readMe << "Rescaled Viscosity =" << VISCOSITY/sqrt(initE*(numOfXGrid-1)*(numOfYGrid-1));
	readMe.close();


	/*======================================================
	Initializing Adams array.
	======================================================*/
	Adams_v = new double**[3];
	Adams_w = new double**[3];
	for(int i = 0;i < 3; i++){
		Adams_v[i] = new double*[numOfXGrid];
		Adams_w[i] = new double*[numOfXGrid];
		for(int j = 0; j < numOfXGrid; j++){
			Adams_v[i][j] = new double[numOfYGrid];
			Adams_w[i][j] = new double[numOfYGrid];
			for(int k = 0; k < numOfYGrid; k++){
				Adams_v[i][j][k] = 0;
				Adams_w[i][j][k] = 0;
			}
		}
	}

	/*==================================================*/
	return;
	
}

Solver_FFTW::~Solver_FFTW(){
}


//function to calculate the first order derivatives
void Solver_FFTW::firstDerivative(){

	/*=============================================
	         calculate dv/dx
	=============================================*/
	for(int i = 0; i < numOfYGrid/2 + 1; i++){
		for(int j = 0; j < numOfXGrid; j++){
			if(j < numOfXGrid/2){
				temp[j][i][0] = -2*PI*j/numOfXGrid*V[j][i][1];
				temp[j][i][1] = 2*PI*j/numOfXGrid*V[j][i][0];
			}
			if(j > numOfXGrid/2){
				temp[j][i][0] = -2*PI*(j-numOfXGrid)/numOfXGrid*V[j][i][1];
				temp[j][i][1] = 2*PI*(j-numOfXGrid)/numOfXGrid*V[j][i][0];
			}
		}
		temp[numOfXGrid/2][i][0] = 0;
		temp[numOfXGrid/2][i][1] = 0;
	}

	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			firstD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			firstD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_firstD);
	for(int i = 0;i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			v_x[i][j] = firstD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}

	/*=================================================
	      calculate dv/dy
	================================================*/
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2; j++){
			temp[i][j][0] = -2*PI*j/numOfYGrid*V[i][j][1];
			temp[i][j][1] = 2*PI*j/numOfYGrid*V[i][j][0];
		}
		temp[i][numOfYGrid/2][0] = 0;
		temp[i][numOfYGrid/2][1] = 0;
	}
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			firstD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			firstD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_firstD);
	for(int i = 0;i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			v_y[i][j] = firstD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}
	/*=================================================
	      calculate dw/dx
	================================================*/
	for(int i = 0; i < numOfYGrid/2 + 1; i++){
		for(int j = 0; j < numOfXGrid; j++){
			if(j < numOfXGrid/2){
				temp[j][i][0] = -2*PI*j/numOfXGrid*W[j][i][1];
				temp[j][i][1] = 2*PI*j/numOfXGrid*W[j][i][0];
			}
			if(j > numOfXGrid/2){
				temp[j][i][0] = -2*PI*(j-numOfXGrid)/numOfXGrid*W[j][i][1];
				temp[j][i][1] = 2*PI*(j-numOfXGrid)/numOfXGrid*W[j][i][0];
			}
		}
		temp[numOfXGrid/2][i][0] = 0;
		temp[numOfXGrid/2][i][1] = 0;
	}

	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			firstD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			firstD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_firstD);
	for(int i = 0;i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			w_x[i][j] = firstD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}
	/*=================================================
	      calculate dw/dy
	==================================================*/
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2; j++){
			temp[i][j][0] = -2*PI*j/numOfYGrid*W[i][j][1];
			temp[i][j][1] = 2*PI*j/numOfYGrid*W[i][j][0];
		}
		temp[i][numOfYGrid/2][0] = 0;
		temp[i][numOfYGrid/2][1] = 0;
	}
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			firstD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			firstD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_firstD);
	for(int i = 0;i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			w_y[i][j] =firstD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}

	/*===============================================*/
	return;
}


//function to calculate the second order derivatives
void Solver_FFTW::secondDerivative(){


	/*=====================================================
	    calculate d^2v/dx^2
	=====================================================*/
	for(int i = 0; i < numOfYGrid/2 + 1;i++){
		for(int j = 0; j < numOfXGrid; j++){
			if(j <= numOfXGrid/2){
				temp[j][i][0] = -(2*PI*j/numOfXGrid)*(2*PI*j/numOfXGrid)*V[j][i][0];
				temp[j][i][1] = -(2*PI*j/numOfXGrid)*(2*PI*j/numOfXGrid)*V[j][i][1];
			}
			if(j > numOfXGrid/2){
				temp[j][i][0] = -(2*PI*(j-numOfXGrid)/numOfXGrid)*(2*PI*(j-numOfXGrid)/numOfXGrid)*V[j][i][0];
				temp[j][i][1] = -(2*PI*(j-numOfXGrid)/numOfXGrid)*(2*PI*(j-numOfXGrid)/numOfXGrid)*V[j][i][1];
			}
		}
	}
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			secondD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			secondD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_secondD);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			v_x_x[i][j] = secondD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}

	/*=====================================================
	    calculate d^2v/dy^2
	=====================================================*/
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2 + 1; j++){
			temp[i][j][0] =  -(2*PI*j/numOfYGrid)*(2*PI*j/numOfYGrid)*V[i][j][0];
			temp[i][j][1] =  -(2*PI*j/numOfYGrid)*(2*PI*j/numOfYGrid)*V[i][j][1];
		}
	}
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			secondD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			secondD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_secondD);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			v_y_y[i][j] = secondD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}

	/*=====================================================
	    calculate d^2w/dx^2
	=====================================================*/
	for(int i = 0; i < numOfYGrid/2 + 1;i++){
		for(int j = 0; j < numOfXGrid; j++){
			if(j <= numOfXGrid/2){
				temp[j][i][0] = -(2*PI*j/numOfXGrid)*(2*PI*j/numOfXGrid)*W[j][i][0];
				temp[j][i][1] = -(2*PI*j/numOfXGrid)*(2*PI*j/numOfXGrid)*W[j][i][1];
			}
			if(j > numOfXGrid/2){
				temp[j][i][0] = -(2*PI*(j-numOfXGrid)/numOfXGrid)*(2*PI*(j-numOfXGrid)/numOfXGrid)*W[j][i][0];
				temp[j][i][1] = -(2*PI*(j-numOfXGrid)/numOfXGrid)*(2*PI*(j-numOfXGrid)/numOfXGrid)*W[j][i][1];
			}
		}
	}
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			secondD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			secondD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_secondD);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			w_x_x[i][j] = secondD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}


	/*=====================================================
	    calculate d^2w/dy^2
	=====================================================*/
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2 + 1; j++){
			temp[i][j][0] =  -(2*PI*j/numOfYGrid)*(2*PI*j/numOfYGrid)*W[i][j][0];
			temp[i][j][1] =  -(2*PI*j/numOfYGrid)*(2*PI*j/numOfYGrid)*W[i][j][1];
		}
	}
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			secondD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			secondD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	fftw_execute(plan_secondD);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			w_y_y[i][j] = secondD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}

	/*===================================================*/
	return;
}


//function of Adams Bashforth method
void Solver_FFTW::adamsMethod(int t){
		//when t=0, Euler's method is used
		if(t == 0){
			for(int i = 0; i < numOfXGrid; i++){
				for(int j = 0; j < numOfYGrid; j++){
					Adams_v[2][i][j] = -(v[i][j]*v_x[i][j] + w[i][j]*v_y[i][j]) + VISCOSITY*(v_x_x[i][j] + v_y_y[i][j]) + externalFx[i][j];
					Adams_w[2][i][j] = -(v[i][j]*w_x[i][j] + w[i][j]*w_y[i][j]) + VISCOSITY*(w_x_x[i][j] + w_y_y[i][j]) + externalFy[i][j];
					v[i][j] = v[i][j] + TIME_STEP*Adams_v[2][i][j];
					w[i][j] = w[i][j] + TIME_STEP*Adams_w[2][i][j];
				}
			}
		}

		//when t=1, second order Adams method is used.
		if(t == 1){
			for(int i = 0; i < numOfXGrid; i++){
				for(int j = 0; j < numOfYGrid; j++){
					Adams_v[1][i][j] = Adams_v[2][i][j];
					Adams_w[1][i][j] = Adams_w[2][i][j];
					Adams_v[2][i][j] = -(v[i][j]*v_x[i][j] + w[i][j]*v_y[i][j]) + VISCOSITY*(v_x_x[i][j] + v_y_y[i][j]) + externalFx[i][j];
					Adams_w[2][i][j] = -(v[i][j]*w_x[i][j] + w[i][j]*w_y[i][j]) + VISCOSITY*(w_x_x[i][j] + w_y_y[i][j]) + externalFy[i][j];
					v[i][j] = v[i][j] + TIME_STEP*(1.5*Adams_v[2][i][j] - 0.5*Adams_v[1][i][j]);
					w[i][j] = w[i][j] + TIME_STEP*(1.5*Adams_w[2][i][j] - 0.5*Adams_w[1][i][j]);
				}
			}
		}

		//when t>=2, third order Adams methods is used.
		if(t >= 2){
			for(int i = 0; i < numOfXGrid; i++){
				for(int j = 0; j < numOfYGrid; j++){
					Adams_v[0][i][j] = Adams_v[1][i][j];
					Adams_w[0][i][j] = Adams_w[1][i][j];
					Adams_v[1][i][j] = Adams_v[2][i][j];
					Adams_w[1][i][j] = Adams_w[2][i][j];
					Adams_v[2][i][j] = -(v[i][j]*v_x[i][j] + w[i][j]*v_y[i][j]) + VISCOSITY*(v_x_x[i][j] + v_y_y[i][j]) + externalFx[i][j];
					Adams_w[2][i][j] = -(v[i][j]*w_x[i][j] + w[i][j]*w_y[i][j]) + VISCOSITY*(w_x_x[i][j] + w_y_y[i][j]) + externalFy[i][j];
					v[i][j] = v[i][j] + TIME_STEP*(23.0/12.0*Adams_v[2][i][j] - 4.0/3.0*Adams_v[1][i][j] + 5.0/12.0*Adams_v[0][i][j]);
					w[i][j] = w[i][j] + TIME_STEP*(23.0/12.0*Adams_w[2][i][j] - 4.0/3.0*Adams_w[1][i][j] + 5.0/12.0*Adams_w[0][i][j]);

				}
			}
		}

		return;
}

// A function to calculate the energy
double Solver_FFTW::calculateE(){
	double E = 0;
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			E += 0.5*(v[i][j]*v[i][j] + w[i][j]*w[i][j]);
		}
	}
	return E;

}


// A function to sample the random force
void Solver_FFTW::samplingForce(){
	for(int i = 0;i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			externalFx[i][j] = 0;
			externalFy[i][j] = 0;
		}
	}
	return;

}

// A function to solver Burgers equation
void Solver_FFTW::burgersSolver_FFTW(){
	/*===============================================
		get V at t=0
	===============================================*/
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			temp_Velocity[i*numOfYGrid+j] = v[i][j];
		}
	}
	fftw_execute(plan_r2c);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2 + 1; j++){
			V[i][j][0] = temp_U[i*(numOfYGrid/2+1) + j][0];
			V[i][j][1] = temp_U[i*(numOfYGrid/2+1) + j][1];
		}
	}
	fftw_execute(plan_c2r);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			v[i][j] = temp_Velocity[i*numOfYGrid + j]/(numOfXGrid*numOfYGrid);
		}
	}
	/*===============================================
		get W at t=0
	===============================================*/
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			temp_Velocity[i*numOfYGrid+j] = w[i][j];
		}
	}
	fftw_execute(plan_r2c);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2 + 1; j++){
			W[i][j][0] = temp_U[i*(numOfYGrid/2+1)+j][0];
			W[i][j][1] = temp_U[i*(numOfYGrid/2+1)+j][1];
		}
	}
	fftw_execute(plan_c2r);
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			w[i][j] = temp_Velocity[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}

	Output *out = new Output(numOfXGrid,numOfYGrid,v,w,0,initE,_VELOCITY);
	/*==============================================
	Time step iteration
	==============================================*/

	for(long int t = 0; t < TIME_N;t++){

		//First step, get V and W
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				temp_Velocity[i*numOfYGrid+j] = v[i][j];
			}
		}
		fftw_execute(plan_r2c);
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid/2 + 1; j++){
				V[i][j][0] = temp_U[i*(numOfYGrid/2+1)+j][0];
				V[i][j][1] = temp_U[i*(numOfYGrid/2+1)+j][1];
			}
		}
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				temp_Velocity[i*numOfYGrid+j] = w[i][j];
			}
		}
		fftw_execute(plan_r2c);
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid/2 + 1; j++){
				W[i][j][0] = temp_U[i*(numOfYGrid/2+1)+j][0];
				W[i][j][1] = temp_U[i*(numOfYGrid/2+1)+j][1];
			}
		}

		/*=============================
		Second step: get v_x,v_y,w_x,w_y
		/*=============================*/
		firstDerivative();
		/*
		if(t == 0){
			Output* out_1 = new Output(numOfXGrid,numOfYGrid,v_x,v_y,0,initE,_DERIVATIVEv);
			Output* out_2 = new Output(numOfYGrid,numOfYGrid,w_x,w_y,0,initE,_DERIVATIVEw);
		}
		*/

		/*=============================
		Third step: get v_x_x,v_y_y,w_x_x,w_y_y
		=============================*/
		secondDerivative();
		/*
		if(t == 0){
			Output* out_1 = new Output(numOfXGrid,numOfYGrid,v_x_x,v_y_y,0,initE,_DDERIVATIVEv);
			Output* out_2 = new Output(numOfYGrid,numOfYGrid,w_x_x,w_y_y,0,initE,_DDERIVATIVEw);
		}
		*/
		
		/*=============================
		Forth step: sampling force
		=============================*/
		samplingForce();
		
		/*=============================
		Last step: Adams-Bashforth method
		=============================*/
		adamsMethod(t);
		
		
		//calculate energy in some steps, this energy is not rescaled
		if((t+1)%ENERGY_OUTPUT == 0){
			double E = calculateE();
			energy << log((t+1)*TIME_STEP) + log(sqrt(initE*(numOfYGrid-1)/(numOfXGrid-1))/(numOfXGrid-1)) << "\t";
			energy << log(E) - log(initE*(numOfYGrid-1)/(numOfXGrid-1))<< endl;
		}

		//generate output in some steps, this output is rescaled, see output.cpp
		if((t+1)%GENERATE_OUTPUT == 0){
			Output* out_1 = new Output(numOfXGrid,numOfYGrid,v,w,t+1,initE,_VELOCITY);
			Output* out_2 = new Output(numOfXGrid,numOfYGrid,v_x,v_y,t+1,initE,_DERIVATIVEv);
			Output* out_3 = new Output(numOfYGrid,numOfYGrid,w_x,w_y,t+1,initE,_DERIVATIVEw);
			cout << "t=" << t+1 << "_completed" << endl;	
		}

	}

	//empty the memory
	energy.close();
	fftw_destroy_plan(plan_c2r);
	fftw_destroy_plan(plan_r2c);
	fftw_destroy_plan(plan_firstD);
	fftw_destroy_plan(plan_secondD);
	return;
}