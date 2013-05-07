/****************************
This is the CPP for output data
******************************/

#include <output.h>

Output::Output(int numOfXGrid, int numOfYGrid, double **xVelocity, double **yVelocity, long int time,double energy){
	ofstream outfile;
	stringstream outfileName;
	outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << ".txt";
	outfile.open(outfileName.str().c_str());
	outfile << time*TIME_STEP*sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1))/(numOfXGrid-1) << endl;
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			outfile << (double)i/(numOfXGrid-1) << "\t" << (double)j/(numOfYGrid-1) << "\t";
			outfile << xVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << "\t" << yVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << endl;
		}
	}
	outfile.close();
	return;
}

Output::~Output(){
}