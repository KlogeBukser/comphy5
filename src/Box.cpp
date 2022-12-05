
#include <iostream>
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

Box::Box(int M, int n_time, double xc, double yc, double px, double py, double sig_x, double sig_y) {
	/*
	Initializer for Box class
	*/
	dx = 1.0 / M;							// Spatial step length (x and y)
	x_dim = (M - 2);						// Number of points in x-direction (and y)
	sq_dim = x_dim * x_dim;					// Number of points in each column/row of A and B
	dt = 1.0 / n_time;						// Time step length

	double v0 = 1.0e5;
	make_potential(v0);
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

	double x = 0.0;
	double y = 0.0;

	u = cx_vec(sq_dim);
	for (int i = 0; i < x_dim; i++) {
		x += dx;
		for (int j = 0; j < x_dim; j++) {
			y += dx;
			u[convert_indices(i, j)] = exp(-0.5 * (pow((x - xc) / sig_x, 2) + pow((y - yc) / sig_y, 2)) + 1i * (px * (x - xc) + py * (y - yc)));
		}
		y = 0.0;
	}
	u = normalise(u);
}

void Box::update_state() {
	cx_vec b = B * u;
	u = spsolve(A, b);
}

void Box::make_potential(double v0) {
	V = mat(sq_dim, sq_dim, fill::value(0.));
}
