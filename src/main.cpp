

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <armadillo>

#include "free_functions.h"



#define ERR_MISSING_PARAMS 1

// Standard library
using namespace std;
// Armadillo
using namespace arma;
// Free functions
using namespace free_funcs;

int main(int argc, char* argv[])
{
    // Random seed
    srand(time(NULL));

    if (!eval_argc(argc)) {
        return ERR_MISSING_PARAMS;
    }

    std::string instruction = argv[1];

    if (instruction == "test") {
        
            
        int M = atoi(argv[2]);                  // Number of spatial steps
        int n_time = atoi(argv[3]);

        double h = 1.0 / M;             // Spatial step length
        double dt = 1.0 / n_time;               // Time step length

        int n_points = (M - 2) * (M - 2);       // Number of points in each column/row of A and B
        

        cx_mat A(n_points,n_points, fill::value(0.));
        cx_mat B(n_points, n_points, fill::value(0.));
        cx_mat V(n_points, n_points, fill::value(0.));

        fill_matrices(A, B, V, h, dt, M);

        A.print();
        B.print();


        // Makes filename
        string filename = make_filename(instruction);
        cout << filename << endl;

    }
}