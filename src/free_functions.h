#ifndef FREE_FUNCTIONS_H
#define FREE_FUNCTIONS_H


namespace free_funcs {
	bool eval_argc(int argc);
	void wrong_cmd_input(std::string instruction);
	std::string make_filename(std::string instruction);

	int convert_indices(int i, int j, int n_spatial);
	arma::sp_cx_mat make_CN_matrix(arma::cx_double r, arma::cx_vec b);
	void fill_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::mat V, double h, double dt, int n_spatial);

	void sp_print(arma::sp_cx_mat& A);
	void fill_potential_mat(arma::mat& V,double v0);
	void update_state(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::cx_vec& u);
	arma::cx_vec set_initial_state(int M, double h, double xc, double yc, double px, double py, double sig_x, double sig_y);
};
#endif