#include <iostream>
#include <string>
#include <armadillo>

#include "free_functions.h"

using namespace std;
using namespace arma;

namespace free_funcs {

	bool eval_argc(int argc) {
		/*
		Evaluates if the minimum amount of cmd parameters are included. Writes information about missing cmd parameters if not enough

		Input: (int) number of cmd parameters
		Output: bool
		*/

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
	
	void wrong_cmd_input(string instruction) {
		/* 
		Writes information about missing cmd parameters (In the case where a correct instruction parameter was given)

		Input: (string) instruction
		Output: void
		*/

		// Is called regardless of what is missing
		cout << "To run program with instruction '" << instruction << "', the following parameters : are needed\n"
			"- Run time instruction: " << instruction << "\n"
			"- number of points in x/y direction (int) \n"
			"- Number of timesteps (int) \n";
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

	int convert_indices(int i, int j, int n_spatial) {
		/*
		Computes index k in flattened matrix from indices i,j from non-flat square matrix

		Input: (int) i,j,n_spatial (i,j:matrix indeces. n_spatial: points in on col/row of old matrix)
		Output: (int) k (index in flattened matrix)
		*/
		return i*n_spatial + j;
	}

	sp_cx_mat make_CN_matrix(cx_double r, cx_vec b) {
		/*
		Produces sparse complex matrix B with vector b as the diagonal and double r on some selected sub/super-diagonals

		Input:
		- complex vector b
		- double value r

		Output: sparse complex matrix B
		*/

		sp_cx_mat B = sp_cx_mat(diagmat(b));
		int sub_dim = sqrt(b.size());

		// Declaration for combined index k
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

	void fill_matrices(sp_cx_mat& A, sp_cx_mat& B, mat V, double h, double dt, int M){
		/*
		- Computes r, a, b from matrix V
		- Fills in matrices A, B from these values

		Input: 
		- References to sparse complex matrices A and B (Used for updating wavefunction in time)
		- Real matrix V (potential values)
		- (double) h,dt (distance between points,smallest time-interval)
		- (int) M (number of points in spatial directions x and y)

		Output: void
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