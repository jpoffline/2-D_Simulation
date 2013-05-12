/****************************
This is the CPP for output data

Author: Francis Chen
Date: 04.03.2013
******************************/

#include <output.h>

Output::Output(int numOfXGrid, int numOfYGrid, double **xVelocity, double **yVelocity, long int time,double energy,output_Type fileType){
	ofstream outfile;
	stringstream outfileName;

	//output file name: Nx= _Ny= _t= .txt
	switch(fileType){
	case _VELOCITY:{
			outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << "_ReT=";
			outfileName << time*TIME_STEP*sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1))/(numOfXGrid-1) << "_VELOCITY.txt";
			outfile.open(outfileName.str().c_str());
		/*================================================================================
		code for the output after rescaling
		=================================================================================*/
		//inputting rescaled time.
		
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				outfile << (double)i/(numOfXGrid-1) << "\t" << (double)j/(numOfYGrid-1) << "\t";
				outfile << xVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << "\t" << yVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << endl;
			}
		}
		/*================================================================================
		code for the output before rescaling
		=================================================================================*/
		/*
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				outfile << i << "\t" << j << "\t" << xVelocity[i][j] << "\t" << yVelocity[i][j] << endl;
			}
		}
		*/
		break;
			}

	//store the derivative of v_x and v_y, the rescaling of derivatives are different from velocities
	case _DERIVATIVEv:{
		outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << "_ReT=";
		outfileName << time*TIME_STEP*sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1))/(numOfXGrid-1) << "_DERIVATIVEv.txt";
		outfile.open(outfileName.str().c_str());
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				outfile << (double)i/(numOfXGrid-1) << "\t" << (double)j/(numOfYGrid-1) << "\t";
				outfile << (numOfXGrid-1)*xVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << "\t" << (numOfXGrid-1)*yVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << endl;
			}
		}
		break;

					 }
	//store the derivative of w_x and w_y, the rescaling of derivatives are different from velocities
	 case _DERIVATIVEw:{
		outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << "_ReT=";
		outfileName << time*TIME_STEP*sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1))/(numOfXGrid-1) << "_DERIVATIVEw.txt";
		outfile.open(outfileName.str().c_str());
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				outfile << (double)i/(numOfXGrid-1) << "\t" << (double)j/(numOfYGrid-1) << "\t";
				outfile << (numOfXGrid-1)*xVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << "\t" << (numOfXGrid-1)*yVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << endl;
			}
		}
		break;
					 }

	//store the derivative of v_x_x and v_y_y
	 case _DDERIVATIVEv:{
		outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << "_ReT=";
		outfileName << time*TIME_STEP*sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1))/(numOfXGrid-1) << "_DDERIVATIVEv.txt";
		outfile.open(outfileName.str().c_str());
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				outfile << (double)i/(numOfXGrid-1) << "\t" << (double)j/(numOfYGrid-1) << "\t";
				outfile << (numOfXGrid-1)*(numOfXGrid-1)*xVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << "\t" << (numOfXGrid-1)*(numOfXGrid-1)*yVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << endl;
			}
		}
		break;
						}

	//store the derivative of w_x_x and w_y_y
	 case _DDERIVATIVEw:{
		outfileName << OUTPUT_PATH << "Nx=" << numOfXGrid << "_Ny=" << numOfYGrid << "_t=" << time << "_ReT=";
		outfileName << time*TIME_STEP*sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1))/(numOfXGrid-1) << "_DDERIVATIVEw.txt";
		outfile.open(outfileName.str().c_str());
		for(int i = 0; i < numOfXGrid; i++){
			for(int j = 0; j < numOfYGrid; j++){
				outfile << (double)i/(numOfXGrid-1) << "\t" << (double)j/(numOfYGrid-1) << "\t";
				outfile << (numOfXGrid-1)*(numOfXGrid-1)*xVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << "\t" << (numOfXGrid-1)*(numOfXGrid-1)*yVelocity[i][j]/sqrt(energy*(numOfYGrid-1)/(numOfXGrid-1)) << endl;
			}
		}
		break;
						}
	//default case, program stops.
	default:{
		cout << "output error!" << endl;
		exit(0);
			}
	}

	outfile.close();
	return;
}

Output::~Output(){
}