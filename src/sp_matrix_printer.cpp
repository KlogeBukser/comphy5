#include <iomanip>
#include <iostream>
#include <armadillo>
#include <assert.h>
#include <complex>

void sp_print(arma::sp_cx_mat &A){
    using namespace std;
    using namespace arma;

    sp_cx_mat::const_iterator it     = A.begin();
    sp_cx_mat::const_iterator it_end = A.end();

    char S[A.n_rows][A.n_cols];


    for (int i =0; i < A.n_rows; i++){
        for (int j = 0; j < A.n_cols; j++){
            S[i][j] = ' ';
        }
    }

    int nnz = 0;
    for(; it != it_end; ++it){
        S[it.row()][it.col()] = '.';
        nnz++;
    }


    for (int i =0; i < A.n_rows; i++){
        cout << "| ";
        for (int j = 0; j < A.n_cols; j++){
            cout << S[i][j] << " ";
        }
        cout <<  "|\n";
    }
    cout << "[matrix size: " << A.n_rows << "x" << A.n_cols << "; n_nonzero: " << nnz << ";]\n";
}
