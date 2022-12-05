#pragma once
#ifndef BOX_H
#define BOX_H


class Box {

	int x_dim;
	int sq_dim;
	double dx;
	double dt;

	arma::cx_vec u;

	arma::mat V;
	arma::sp_cx_mat A;
	arma::sp_cx_mat B;

public:
	Box(int M, int n_time, double xc, double yc, double px, double py, double sig_x, double sig_y);

	int convert_indices(int i, int j);
	void make_matrices();

	void make_potential(double v0);
	void print();
	
	void set_initial_state(double xc, double yc, double px, double py, double sig_x, double sig_y);
	void update_state();
};
#endif