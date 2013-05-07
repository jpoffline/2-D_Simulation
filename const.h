/*******************************************************
 This is the constants header file in simulation
 the 2-D model of Burgers equation        
 
 Author: Francis Chen
 Date: 04.03.2013
 *******************************************************/ 

#define VISCOSITY 1000 // This is the viscocity \nu in the equation
#define TIME_N 400000     // This is the total time steps in the simulation
#define TIME_STEP 1e-5    // This is the delta t in simulation
#define PI 3.1415926535897932 
#define INPUT_PATH "D:/Heidelberg/2-D_Executives/input/" //This is the input path
#define OUTPUT_PATH "D:/Heidelberg/2-D_Executives/output/" // This is the input path
#define ENERGY_OUTPUT 10000 // after how many iterations one want to calculate the energy
#define GENERATE_OUTPUT 20000  //after how many iterations one want to generate an output
#define THREADS 2 //how many threads used in the simulation