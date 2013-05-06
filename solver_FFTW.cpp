/********************************
This is the CPP file of FFTW solver
**********************************/

#include "solver_FFTW.h"

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
	temp_Velocity = new double[numOfXGrid*numOfYGrid];
	firstD_u = new double[numOfXGrid*numOfYGrid];
	secondD_u = new double[numOfXGrid*numOfYGrid];
	for(int i = 0; i < numOfXGrid; i++){
		v[i] = new double[numOfYGrid];
		w[i] = new double[numOfYGrid];
		for(int j = 0; j < numOfYGrid; j++){
			v[i][j] = initInput->getXVelocity(i,j);
			w[i][j] = initInput->getYVelocity(i,j);
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
	initializing multiple threads
	==========================================================*/
	if(fftw_init_threads()){
		fftw_plan_with_nthreads(THREADS);
		cout << "Using "<< THREADS << " threads" << endl;
	}
	else {
		cout << "Using multiple threads failed" << endl;
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

	Output *out = new Output(numOfXGrid,numOfYGrid,v,w,0);
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

		//Second step: get v_x,v_y,w_x,w_y
		firstDerivative();

		//Third step: get v_x_x,v_y_y,w_x_x,w_y_y
		secondDerivative();

		//Last step: Adams-Bashforth method
		if(t == 0){
			for(int i = 0; i < numOfXGrid; i++){
				for(int j = 0; j < numOfYGrid; j++){
					Adams_v[2][i][j] = -(v[i][j]*v_x[i][j] + w[i][j]*v_y[i][j]) + VISCOSITY*(v_x_x[i][j] + v_y_y[i][j]);
					Adams_w[2][i][j] = -(v[i][j]*w_x[i][j] + w[i][j]*w_y[i][j]) + VISCOSITY*(w_x_x[i][j] + w_y_y[i][j]);
					v[i][j] = v[i][j] + TIME_STEP*Adams_v[2][i][j];
					w[i][j] = w[i][j] + TIME_STEP*Adams_w[2][i][j];
				}
			}
		}
		if(t == 1){
			for(int i = 0; i < numOfXGrid; i++){
				for(int j = 0; j < numOfYGrid; j++){
					Adams_v[1][i][j] = Adams_v[2][i][j];
					Adams_w[1][i][j] = Adams_w[2][i][j];
					Adams_v[2][i][j] = -(v[i][j]*v_x[i][j] + w[i][j]*v_y[i][j]) + VISCOSITY*(v_x_x[i][j] + v_y_y[i][j]);
					Adams_w[2][i][j] = -(v[i][j]*w_x[i][j] + w[i][j]*w_y[i][j]) + VISCOSITY*(w_x_x[i][j] + w_y_y[i][j]);
					v[i][j] = v[i][j] + TIME_STEP*(1.5*Adams_v[2][i][j] - 0.5*Adams_v[1][i][j]);
					w[i][j] = w[i][j] + TIME_STEP*(1.5*Adams_w[2][i][j] - 0.5*Adams_w[1][i][j]);
				}
			}
		}
		if(t >= 2){
			for(int i = 0; i < numOfXGrid; i++){
				for(int j = 0; j < numOfYGrid; j++){
					Adams_v[0][i][j] = Adams_v[1][i][j];
					Adams_w[0][i][j] = Adams_w[1][i][j];
					Adams_v[1][i][j] = Adams_v[2][i][j];
					Adams_w[1][i][j] = Adams_w[2][i][j];
					Adams_v[2][i][j] = -(v[i][j]*v_x[i][j] + w[i][j]*v_y[i][j]) + VISCOSITY*(v_x_x[i][j] + v_y_y[i][j]);
					Adams_w[2][i][j] = -(v[i][j]*w_x[i][j] + w[i][j]*w_y[i][j]) + VISCOSITY*(w_x_x[i][j] + w_y_y[i][j]);
					v[i][j] = v[i][j] + TIME_STEP*(23.0/12.0*Adams_v[2][i][j] - 4.0/3.0*Adams_v[1][i][j] + 5.0/12.0*Adams_v[0][i][j]);
					w[i][j] = w[i][j] + TIME_STEP*(23.0/12.0*Adams_w[2][i][j] - 4.0/3.0*Adams_w[1][i][j] + 5.0/12.0*Adams_w[0][i][j]);

				}
			}
		}

		//generate some outputs
//		Output* out = new Output(numOfXGrid,numOfYGrid,v,w,t+1);
		
		if((t+1)%ENERGY_OUTPUT == 0){
			double E = calculateE();
			energy << log((t+1)*TIME_STEP) << "\t" << log(E) << endl;
		}
		if((t+1)%GENERATE_OUTPUT == 0){
			Output* out = new Output(numOfXGrid,numOfYGrid,v,w,t+1);
			cout << "t=" << t+1 << "_completed" << endl;	
		}

	}


	energy.close();
	fftw_destroy_plan(plan_c2r);
	fftw_destroy_plan(plan_r2c);
	fftw_destroy_plan(plan_firstD);
	fftw_destroy_plan(plan_secondD);
	return;
}