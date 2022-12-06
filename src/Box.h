#pragma once
#ifndef BOX_H
#define BOX_H


class Box {

	int x_dim;
	int sq_dim;
	double dx;
	double dt;
	double current_time = 0.0;

	arma::cx_vec u;
	arma::vec x;
	arma::vec y;

	arma::mat V;
	arma::sp_cx_mat A;
	arma::sp_cx_mat B;

public:
	Box(int n_slits, double h, double timestep, double xc, double sig_x, double px, double yc, double sig_y, double py, double V0);
	
	int convert_indices(int i, int j);
	void make_pos_vectors();
	void make_matrices();
	void make_potential(int n_slits, double v0);
	void make_potential_double(double v0);

	void print();
	std::string get_string();
	double total_probability();

	void set_initial_state(double xc, double yc, double px, double py, double sig_x, double sig_y);
	void update_state();
	void evolve_probability(double time,std::string filename);
};
#endif