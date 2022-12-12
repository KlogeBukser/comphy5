
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
	x_dim = int(1.0 / dx) - 1;				// (M - 2) = Number of points in x-direction (and y)
	sq_dim = x_dim * x_dim;					// Number of points in each column/row of A and B
	dt = timestep;							// Time step length

	

	make_pos_vectors();
	make_potential(n_slits,V0);
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
	
	superlu_opts opts;

	opts.allow_ugly = true;
	opts.symmetric = true;

	spsolve(u, A, b, "superlu", opts);
	current_time += dt;
}

void Box::update_state(double time) {
	// Updates the current time until the desired time is reached
	while (current_time <= time - 1e-10) {
		update_state();
	}
}


string Box::get_string(int i, int j) {
	int k = convert_indices(i, j);

	return to_string(norm(u[k])) + " " + to_string(real(u[k])) + " " + to_string(imag(u[k]));
}

double Box::total_probability() {
	
	double probability = 0.0;
	for (int k = 0; k < sq_dim; k++) {
		probability += norm(u[k]);
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
	
	V = vec(sq_dim);
	if (n_slits == 2) {
		make_potential_double(v0);
	}
}

void Box::make_potential_double(double v0) {
	
	double thickness = 0.02;

	for (int i = 0; i < x_dim; i++) {
		if (abs(x[i] - 0.5) > thickness / 2.0) {
			continue;
		}
		for (int j = 0; j < x_dim; j++) {

			// Skips the current point if it is in the slit opening
			if (!(abs(0.45 - y[j]) < 0.025 || abs(0.55 - y[j]) < 0.025)) {
				V[convert_indices(i, j)] = v0;
			}
		}
	}
	
}

void Box::write_probability(double time, string filename) {

	ofstream outfile(filename);

	// Initial probability (should be 1)
	outfile << current_time << " , " << to_string(total_probability()) << endl;

	while (current_time < time) {
		update_state();
		outfile << current_time << " , " << to_string(total_probability()) << endl;	
		cout << "Time: " << current_time << endl;
	}
}

void Box::write_state(vec times, string filename) {

	ofstream outfile(filename);

	outfile << to_string(x_dim) << endl;

	for (double time : times) {
		update_state(time);
		outfile << to_string(current_time) << " : " << "0 0 " << get_string(0, 0);

		for (int i = 1; i < x_dim; i++) {
			for (int j = 1; j < x_dim; j++) {
				outfile << " , " << to_string(i) << " " << to_string(j) << " " << get_string(i, j);
			}
		}
		outfile << endl;
	}
	outfile.close();

	ofstream outfile2("textfiles/interference.txt");

	for (int i = 0; i < x_dim; i++) {
		if (abs(x[i] - 0.8) < 1e-8) {

			for (int j = 0; j < x_dim; j++) {
				outfile2 << to_string(y[j]) << " " << to_string(norm(u[convert_indices(i, j)])) << endl;
			}

		}
	}
	outfile2.close();

}

void Box::write_state(double time, string filename) {

	ofstream outfile(filename);
	
	outfile << to_string(x_dim) << endl;

	while (current_time < time) {
		outfile << to_string(current_time) << " : " << "0 0 " << get_string(0, 0);

		for (int i = 1; i < x_dim; i++) {
			for (int j = 1; j < x_dim; j++) {
				outfile << " , " << to_string(i) << " " << to_string(j) << " " << get_string(i, j);
			}
		}
		outfile << endl;

		update_state();
	}
	outfile.close();
}

