
#ifndef BOX_H
#define BOX_H


class Box {

private:

private:

	typedef void (Box::* potential_pointer)(int i, double v0);
	potential_pointer potential;

	int x_dim;
	int sq_dim;
	double dx;
	double dt;
	double current_time = 0.0;
	int n_slits;

	arma::cx_vec u;
	arma::vec x;
	arma::vec y;

	arma::vec V;
	arma::sp_cx_mat A;
	arma::sp_cx_mat B;

public:
	Box(int nunber_of_slits, double h, double timestep, double xc, double sig_x, double px, double yc, double sig_y, double py, double V0);
	
	int convert_indices(int i, int j);
	void make_pos_vectors();
	void make_matrices();
	void make_potential(double v0);
	
	void single_potential(int i, double v0);
	void double_potential(int i, double v0);
	void triple_potential(int i, double v0);

	void print();
	double total_probability();

	void set_initial_state(double xc, double yc, double px, double py, double sig_x, double sig_y);
	void update_state();
	void update_state(double time);
	void write_probability(double time,std::string filename);
	void write_state(arma::vec times, std::string filename);
	void write_state(double time, std::string filename);
};
#endif