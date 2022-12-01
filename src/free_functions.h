#ifndef FREE_FUNCTIONS_H
#define FREE_FUNCTIONS_H


namespace free_funcs {
	bool eval_argc(int argc);
	void wrong_cmd_input(std::string instruction);
	std::string make_filename(std::string instruction);

	int convert_indices(int i, int j, int n_spatial);
	arma::cx_mat make_CN_matrix(arma::cx_double r, arma::cx_vec b);
	void fill_matrices(arma::cx_mat& A, arma::cx_mat& B, arma::cx_mat V, double h, double dt, int n_spatial);
};
#endif