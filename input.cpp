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
	infileName << INPUT_PATH << fn;
	infile.open(infileName.str().c_str());
	if(!infile){
		cerr << "Error: unable to open input file: " << infile << endl;
		exit(0);
	}

	infile >> numOfXGrid;
	infile >> numOfYGrid;
	xVelocity = new double*[numOfXGrid];
	yVelocity = new double*[numOfXGrid];
	for(int i = 0; i < numOfXGrid; i++){
		xVelocity[i] = new double[numOfYGrid];
		yVelocity[i] = new double[numOfYGrid];
	}

	int x, y;
	for(int i = 0; i < numOfXGrid*numOfYGrid; i++){
		infile >> x;
		infile >> y;
		infile >> xVelocity[x][y];
		infile >> yVelocity[x][y];
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

double Input::getXVelocity(int i, int j){
	return xVelocity[i][j];
}

double Input::getYVelocity(int i, int j){
	return yVelocity[i][j];
}