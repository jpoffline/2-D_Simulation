/******************************
This is the CPP file of input data
*****************************/

#include "input.h"

Input::Input(){
	ifstream infile;
	stringstream infileName;
	string fn;
	cout << "Input the file name" << endl;
	cin >> fn;
	infileName >> INPUT_PATH >> fn;
	infile.open(infileName.str().c_str());
	if(!infile){
		cerr << "Error: unable to open input file: " << infile << endl;
	}

	infile >> numOfXGrid;
	infile >> numOfYGrid;
	initVelocity = new double*[numOfXGrid];
	for(int i = 0; i < numOfXGrid; i++){
		initVelocity[i] = new double[numOfYGrid];
		for(int j = 0; j < numOfYGrid; j++){
			infile >> initVelocity[i][j];
		}
	}
	infile.close();
	return;
}

Input::~Input(){
}

int Input::getXGridNum(){
	return numOfXGrid;
}

int Input::getYGridNum(){
	return numOfYGrid;
}

double Input::getVelocity(int i, int j){
	return initVelocity[i][j];
}