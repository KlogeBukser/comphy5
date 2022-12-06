

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <armadillo>

#include "Box.h"
#include "free_functions.h"

// Standard library
using namespace std;
// Armadillo
using namespace arma;


Box::Box(int n_slits, double h, double timestep, double xc, double sig_x, double px, double yc, double sig_y, double py, double V0) {
	/*
	Initializer for Box class
	*/
	dx = h;									// Spatial step length (x and y)
	x_dim = int(1.0 / dx) - 2;				// (M - 2) = Number of points in x-direction (and y)
	sq_dim = x_dim * x_dim;					// Number of points in each column/row of A and B
	dt = timestep;							// Time step length


	make_pos_vectors();
	make_potential(n_slits,V0);
	make_matrices();

	set_initial_state(xc, yc, px, py, sig_x, sig_y);
}

int Box::convert_indices(int i, int j) {
	return i * x_dim + j;
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

	for (int i = 0; i < x_dim; i++) {
		
		for (int j = 0; j < x_dim; j++) {
			int k = convert_indices(i, j);
			b(k) = 1. - 4. * r - r_h_sq * V(i, j);
		}
	}

	B = sp_cx_mat(diagmat(b));

	// Fills the values on sub- and superdiagonal
	for (int i = 0; i < x_dim; i++) {
		for (int j = 0; j < x_dim - 1; j++) {

			int k = convert_indices(i, j);
			B(k, k + 1) = r;
			B(k + 1, k) = r;
		}
	}

	// Fills the values on the diagonals 'x_dim' elements from the main diagonal.
	for (int k = 0; k < x_dim * (x_dim - 1); k++) {
		B(k, k + x_dim) = r;
		B(k + x_dim, k) = r;
	}

	A = 2.0*A.eye(size(B)) - B;
}



void Box::print() {

	cout << "-------- A --------\n\n";
	free_funcs::sp_print(A);
	cout << "-------- B --------\n\n";
	free_funcs::sp_print(B);
}


void Box::set_initial_state(double xc, double yc, double px, double py, double sig_x, double sig_y) {

	u = cx_vec(sq_dim);
	for (int i = 0; i < x_dim; i++) {
		for (int j = 0; j < x_dim; j++) {
			u[convert_indices(i, j)] = exp(-0.5 * (pow((x[i] - xc) / sig_x, 2) + pow((y[j] - yc) / sig_y, 2)) + 1i * (px * (x[i] - xc) + py * (y[j] - yc)));
		}
	}
	u = normalise(u);
}

void Box::update_state() {
	cx_vec b = B * u;
	u = spsolve(A, b);
	current_time += dt;
}

string Box::get_string() {
	string info = to_string(real(u[0])) + " " + to_string(imag(u[0]));
	for (int i = 1; i < sq_dim; i++) {
		info += " , " + to_string(real(u[i])) + " " + to_string(imag(u[i]));
	}
	return info;
}

double Box::total_probability() {
	
	double probability = 0.0;
	for (int i = 0; i < sq_dim; i++) {
		probability += norm(u[i]);
	}
	return probability;
}



void Box::make_pos_vectors() {

	x = vec(x_dim);
	y = vec(x_dim);
	for (int i = 0; i < x_dim; i++) {
		x[i] = dx * (i + 1);
		y[i] = dx * (i + 1);
	}
}

void Box::make_potential(int n_slits, double v0) {
	
	V = mat(x_dim, x_dim, fill::value(0.));
	if (n_slits == 2) {
		make_potential_double(v0);
	}
}

void Box::make_potential_double(double v0) {
	
	double x_mid = 0.5;
	double thickness = 0.02;
	double mid_wall_length = 0.05;
	double aperture = 0.05;

	// This snippet finds which indices in x-direction are inside the wall
	int width_index = 0;
	while (x[width_index] < thickness) {
		width_index++;
	}
	width_index--;

	// The outer loop only loops over the x-indices corresponding to the wall
	for (int i = (x_dim - width_index) / 2; i < (x_dim + width_index) / 2; i++) {
		for (int j = 0; j < x_dim; j++) {

			// Skips the current point if it is in the slit opening
			if (abs(0.45 - y[j]) < 0.025 || abs(0.55 - y[j]) < 0.025) {
				continue;
			}

			V[i, j] = v0;
		}
	}
}

void Box::evolve_probability(double time, string filename) {

	ofstream outfile(filename);

	while (current_time < time) {
		outfile << current_time << " , " << to_string(total_probability()) << endl;
		update_state();
		cout << "Time: " << current_time << endl;
	}
}
