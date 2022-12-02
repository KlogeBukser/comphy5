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
};
#endif