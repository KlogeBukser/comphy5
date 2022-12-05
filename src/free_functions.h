#ifndef FREE_FUNCTIONS_H
#define FREE_FUNCTIONS_H


namespace free_funcs {
	bool eval_argc(int argc);
	void wrong_cmd_input(std::string instruction);
	std::string make_filename(std::string instruction);
	void sp_print(arma::sp_cx_mat& A);
};
#endif