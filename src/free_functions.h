#ifndef FREE_FUNCTIONS_H
#define FREE_FUNCTIONS_H


namespace free_funcs {
	void print_minimum_req();
	bool eval_argc(int argc);
	bool eval_argc(int argc, std::string instruction);
	void wrong_cmd_input(std::string instruction);
	std::string make_filename(std::string instruction);
	std::string make_filename(std::string instruction, int n_slits);
	void sp_print(arma::sp_cx_mat& A);
};
#endif