
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <armadillo>

#include "headers/Box.h"
#include "headers/free_functions.h"

// Standard library
using namespace std;
// Armadillo
using namespace arma;


Box::Box(int nunber_of_slits, double h, double timestep, double xc, double sig_x, double px, double yc, double sig_y, double py, double V0) {
	/*
	Initializer for Box class
	*/
	dx = h;									// Spatial step length (x and y)
	x_dim = int(1.0 / dx) - 1;				// (M - 2) = Number of points in x-direction (and y)
	sq_dim = x_dim * x_dim;					// Number of points in each column/row of A and B
	dt = timestep;							// Time step length
	n_slits = nunber_of_slits;
	

	make_pos_vectors();
	make_potential(V0);
	make_matrices();
	set_initial_state(xc, yc, px, py, sig_x, sig_y);
}

int Box::convert_indices(int i, int j) {
	return j * x_dim + i;
}

void Box::make_matrices() {
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
	cx_double r = r_h_sq / (dx * dx);

	cx_vec b(sq_dim);


	for (int k = 0; k < sq_dim; k++) {
		b(k) = 1. - 4. * r - r_h_sq * V[k];
	}
	B = sp_cx_mat(diagmat(b));

	// Fills the values on sub- and superdiagonal
	for (int k = 0; k < sq_dim - 1; k++) {
		B(k, k + 1) = r;
		B(k + 1, k) = r;
	}

	// Fills the values on the outer diagonals
	for (int k = 0; k < x_dim * (x_dim - 1); k++) {
		B(k, k + x_dim) = r;
		B(k + x_dim, k) = r;
	}

	A = 2.0*A.eye(size(B)) - B;
}



void Box::print() {
	/*
	Prints the sparse matrix, with dots representing the non-zero elements. This is only informative for small matrices.
	*/

	cout << "-------- A --------\n\n";
	free_funcs::sp_print(A);
	cout << "-------- B --------\n\n";
	free_funcs::sp_print(B);
}


void Box::set_initial_state(double xc, double yc, double px, double py, double sig_x, double sig_y) {
	/*
	Sets the initial quantum state based on input parameters for center position, width and momentum in x- and y-direction (input parameters)
	*/

	u = cx_vec(sq_dim);
	for (int i = 0; i < x_dim; i++) {
		for (int j = 0; j < x_dim; j++) {
			u[convert_indices(i, j)] = exp(-0.5 * (pow((x[i] - xc) / sig_x, 2) + pow((y[j] - yc) / sig_y, 2)) + 1i * (px * (x[i] - xc) + py * (y[j] - yc)));
		}
	}
	u = normalise(u);
}

void Box::update_state() {
	/*
	Updates the system state, and increases the time by one step.
	*/
	cx_vec b = B * u;
	
	spsolve(u, A, b);
	current_time += dt;
}

void Box::update_state(double time) {
	/*
	Updates the current time until the desired time is reached, with a small toleranse
	*/ 
	while (current_time < time - 1e-12) {
		update_state();
	}
}


double Box::total_probability() {
	/*
	Sums over the probability density for all points in the box.
	*/
	
	double probability = -1.0;
	for (int k = 0; k < sq_dim; k++) {
		probability +=  norm(u[k]);
	}
	return probability;
}


void Box::make_pos_vectors() {
	/*
	Vectors of the possible position values in x- and y- direction.
	*/

	x = vec(x_dim);
	y = vec(x_dim);
	for (int i = 0; i < x_dim; i++) {
		x[i] = dx * (i + 1);
		y[i] = dx * (i + 1);
	}
}

void Box::make_potential( double v0) {
	/*
	Checks how many slits to use for the potential, and assigns a function pointer to the function for constructing that wall potential.
	Finds the x-indices corresponding to the middle wall, and calls the function-pointer, which loops over the y-values, and fills the potential.
	*/
	V = vec(sq_dim);
	
	if (n_slits == 1) {
		potential = &Box::single_potential;
	}

	else if (n_slits == 2) {
		potential = &Box::double_potential;
	}

	else if (n_slits == 3) {
		potential = &Box::triple_potential;
	}
	else {
		return;
	}

	double thickness = 0.02;
	for (int i = 0; i < x_dim; i++) {
		if (abs(x[i] - 0.5) > thickness / 2.0) {
			continue;
		}

		(this->*potential)(i, v0);
	}
}

void Box::single_potential(int i, double v0) {
	/*
	Potential for a wall with 1 slit.
	Is called from Box::make_potential.
	Is symmetric over y = 0.5, with slit centered at y = 0.5 with 0.05 aperture
	*/

	for (int j = 0; j < x_dim; j++) {

		// Skips the current point if it is in the slit opening
		if (!(abs(0.5 - y[j]) < 0.025)) {
			V[convert_indices(i, j)] = v0;
		}
	}
}

void Box::double_potential(int i, double v0) {
	/*
	Potential for a wall with 2 slits.
	Is called from Box::make_potential.
	Is symmetric over y = 0.5, with slits centered at y = 0.45,0.55, with 0.05 aperture
	*/
	
	for (int j = 0; j < x_dim; j++) {

		// Skips the current point if it is in the slit opening
		if (!(abs(0.45 - y[j]) < 0.025 || abs(0.55 - y[j]) < 0.025)) {
			V[convert_indices(i, j)] = v0;
		}
	}	
}

void Box::triple_potential(int i, double v0) {
	/*
	Potential for a wall with 3 slita.
	Is called from Box::make_potential.
	Is symmetric over y = 0.5, with slits centered at y = 0.4,0.5,0.6 with 0.05 aperture
	*/

	for (int j = 0; j < x_dim; j++) {

		// Skips the current point if it is in the slit opening
		if (!(abs(0.4 - y[j]) < 0.025 || abs(0.5 - y[j]) < 0.025 || abs(0.6 - y[j]) < 0.025)) {
			V[convert_indices(i, j)] = v0;
		}
	}
}

void Box::write_probability(double time, string filename) {
	/*
	Writes the current time and total probability to file and updates the system state for every time-state
	*/

	ofstream outfile(filename);

	// Initial probability (should be 1)
	outfile << current_time << " , " << total_probability() << endl;

	while (current_time < time) {
		update_state();
		outfile << current_time << " , " << total_probability() << endl;
	}
	outfile.close();
}

void Box::write_state(vec times, string filename) {
	/*
	Writes the indices i,j as well as the probability density, real part of density, and imaginary part of density
	for every point to file for specific time values. The system is updated until these time-values are reached.
	*/

	ofstream outfile(filename);

	outfile << x_dim << endl;

	for (double time : times) {
		update_state(time);
		outfile << current_time << " : " << "0 0 " << norm(u[0]) << " " << real(u[0]) << " " << imag(u[0]);

		for (int i = 1; i < x_dim; i++) {
			for (int j = 1; j < x_dim; j++) {
				int k = convert_indices(i, j);
				outfile << " , " << i << " " << j << " " << norm(u[k]) << " " << real(u[k]) << " " << imag(u[k]);
			}
		}
		outfile << endl;
	}
	outfile.close();

	ofstream outfile2("textfiles/interference" + to_string(n_slits) + ".txt");

	outfile2 << 0.8 << " " << current_time << endl;

	for (int i = 0; i < x_dim; i++) {
		if (abs(x[i] - 0.8) < 1e-8) {

			for (int j = 0; j < x_dim; j++) {
				outfile2 << y[j] << " " << norm(u[convert_indices(i, j)]) << endl;
			}

		}
	}
	outfile2.close();

}

void Box::write_state(double time, string filename) {
	/*
	Writes the indices i,j as well as the probability density, real part of density, and imaginary part of density
	for every point to file for every time step. The system is updated until the final time (input parameter) has been reached.
	*/

	ofstream outfile(filename);
	
	outfile << x_dim << endl;

	while (current_time < time) {
		outfile << current_time << " : " << "0 0 " << norm(u[0]) << " " << real(u[0]) << " " << imag(u[0]);

		for (int i = 1; i < x_dim; i++) {
			for (int j = 1; j < x_dim; j++) {
				int k = convert_indices(i, j);
				outfile << " , " << i << " " << j << " " << norm(u[k]) << " " << real(u[k]) << " " << imag(u[k]);
			}
		}
		outfile << endl;

		update_state();
	}
	outfile.close();
}

