#include <iostream>
#include <string>
#include <armadillo>

#include "free_functions.h"

using namespace std;
using namespace arma;

namespace free_funcs {

	bool eval_argc(int argc) {

		if (argc < 3) {
			// These parameters are needed for all program runs
			std::cout << "\nMissing input parameters! (" << argc - 1 << " parameters was included.) \n"
				"All calls of this program should include the following parameters:\n"
				"- Run time instruction: ('simple' / 'burn' / 'distribution' / 'critical' )\n"
				"- number of points in x/y direction (int) \n"
				"- Number of timesteps (int) \n";
			return false;
		}
		return true;
	}
	
	void wrong_cmd_input(std::string instruction) {
		/* Writes error message if the number of arguments is not correct for the given instructions */

		// Is called regardless of what is missing
		std::cout << "To run program with instruction '" << instruction << "', the following parameters : are needed\n"
			"- Run time instruction: " << instruction << "\n"
			"- number of points in x/y direction (int) \n"
			"- Number of timesteps (int) \n";
	}

	std::string make_filename(std::string instruction) {

		std::string filename = "textfiles/" + instruction + ".txt";
		
		return filename;
	}

	int convert_indices(int i, int j, int n_spatial) {
		return i*n_spatial + j;
	}

	arma::cx_mat make_CN_matrix(arma::cx_double r, arma::cx_vec b) {

		cx_mat B = diagmat(b);
		int sub_dim = sqrt(b.size());


		int k;

		// Fills the values on sub- and superdiagonal
		for (int i = 0; i < sub_dim; i++) {
			for (int j = 0; j < sub_dim - 1; j++) {
				
				k = convert_indices(i, j, sub_dim);
				B(k, k + 1) = r;
				B(k + 1, k) = r;
			}
		}
		
		// Fills the values on the diagonals 'sub_dim' elements from the main diagonal.
		for (int k = 0; k < sub_dim * (sub_dim - 1); k++) {
			B(k, k + sub_dim) = r;
			B(k + sub_dim, k) = r;
		}

		return B;
	}

	void fill_matrices(cx_mat& A, cx_mat& B, cx_mat V, double h, double dt, int M){
		/*
		- Computes r, a, b from matrix V
		- Fills in matrices A, B from these values
		*/
		cx_double r_h_sq = 1i * dt / 2.;                  // r*h^2
		cx_double r = r_h_sq / (h * h);
		int n_points = (M - 2) * (M - 2);

		cx_vec a(n_points);
		cx_vec b(n_points);
		int k;
		for (int i = 0; i < M - 2; i++) {
			for (int j = 0; j < M - 2; j++) {
				k = convert_indices(i, j, M-2);
				a(k) = 1. + 4. * r + r_h_sq * V(i, j);
				b(k) = 2. - a(k);
			}
		}

		A = make_CN_matrix(-r, a);
		B = make_CN_matrix(r, b);
	}



}