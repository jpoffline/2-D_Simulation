/****************************
This is the CPP for output data
******************************/

#include <output.h>

Output::Output(int numOfXGrid, int numOfYGrid, double **xVelocity, double **yVelocity, long int time){
	ofstream outfile;
	stringstream outfileName;
	outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << ".txt";
	outfile.open(outfileName.str().c_str());
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			outfile << i << "\t" << j << "\t" << xVelocity[i][j] << "\t" << yVelocity[i][j] << endl;
		}
	}
	outfile.close();
	return;
}

Output::~Output(){
}