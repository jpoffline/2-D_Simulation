/****************************
This is the CPP for output data
******************************/

#include <output.h>

Output::Output(int numOfXGrid, int numOfYGrid, double **tempVelocity, int time){
	ofstream outfile;
	stringstream outfileName;
	outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << ".txt";
	outfile.open(outfileName.str().c_str());
	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			outfile << i << "\t" << j << "\t" << tempVelocity[i][j] << endl;
		}
	}
	outfile.close();
	return;
}

Output::~Output(){
}