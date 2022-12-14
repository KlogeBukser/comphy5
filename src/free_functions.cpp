#include <iostream>
#include <string>
#include <armadillo>

#include "headers/free_functions.h"

using namespace std;
using namespace arma;

namespace free_funcs {

	void print_minimum_req() {
		cout << 
			"All calls of this program should include the following parameters:\n"
			"- Run time instruction: ('probability' / 'heat' / 'evo')\n"
			"- number of slits (int: 0,1,2,3) \n"
			"- step length 'h' in x/y direction (double) \n"
			"- timestep 'dt' (double) \n"
			"- center of wavepacket 'xc' in x-direction (double) \n"
			"- width of wavepacket 'sig_x' in x-direction (double) \n"
			"- momentum of wavepacket 'px' in x-direction (double) \n"
			"- center of wavepacket 'yc' in y-direction (double) \n"
			"- width of wavepacket 'sig_y' in y-direction (double) \n"
			"- momentum of wavepacket 'py' in y-direction (double) \n"
			"- Total time 'T' (double) \n";
	}


	bool eval_argc(int argc) {
		/*
		Evaluates if the minimum amount of cmd parameters are included. Writes information about missing cmd parameters if not enough

		Input: (int) number of cmd parameters
		Output: bool
		*/

		if (argc < 12) {
			// (argc - 1) is used here as the first parameter 'main.out' is not a paramater in the same sense as the others 
			cout << "\nMissing input parameters! (" << argc - 1 << " parameters was included) \n";
			print_minimum_req();
			
			
			return false;
		}
		return true;
	}

	
	void wrong_cmd_input(string instruction) {
		/* 
		Writes information about missing cmd parameters (In the case where a correct instruction parameter was given)

		Input: (string) instruction
		Output: void
		*/

		// Is called regardless of what is missing
		cout << "To run program with instruction '" << instruction << "', the following parameters : are needed:\n";
		print_minimum_req();

	}
	string make_filename(string instruction, int n_slits) {
		return make_filename(instruction + "_" + to_string(n_slits));
	}

	string make_filename(string instruction) {
		/*
		Makes filename for a txt file

		Input: (string) instruction
		Output: (string) filename
		*/

		string filename = "textfiles/" + instruction + ".txt";
		
		return filename;
	}

	void sp_print(sp_cx_mat& A) {

		/*
		From lab assistant Håkon Olav Torvik, via Thomas Bergheim (with slight modifications).
		
		Prints the structure of a complex sparse matrix

		Input: Reference to sparse complex matrix A
		Output: void
		*/

		sp_cx_mat::const_iterator it = A.begin();
		sp_cx_mat::const_iterator it_end = A.end();

		char S[A.n_rows][A.n_cols];

		for (int i = 0; i < A.n_rows; i++) {
			for (int j = 0; j < A.n_cols; j++) {
				S[i][j] = ' ';
			}
		}

		int nnz = 0;
		for (; it != it_end; ++it) {
			S[it.row()][it.col()] = '.';
			nnz++;
		}


		for (int i = 0; i < A.n_rows; i++) {
			cout << "| ";
			for (int j = 0; j < A.n_cols; j++) {
				cout << S[i][j] << " ";
			}
			cout << "|\n";
		}
		cout << "[matrix size: " << A.n_rows << "x" << A.n_cols << "; n_nonzero: " << nnz << ";]\n";
	}
}