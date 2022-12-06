

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <armadillo>

#include "free_functions.h"
#include "Box.h"



constexpr auto ERR_MISSING_PARAMS = 1;

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
    int n_slits = atoi(argv[2]);
    double barrier_potential = 0.0;

    // Takes the cmd_input for potential if there is a barrier in the middle. 
    if (n_slits > 0) {
        if (!eval_argc(argc,"barrier")) {
            return ERR_MISSING_PARAMS;
        }
        barrier_potential = atof(argv[12]);
    }

    double h = atof(argv[3]);
    double dt = atof(argv[4]);
    double time = atof(argv[5]);

    double xc = atof(argv[6]);
    double sig_x = atof(argv[7]);
    double px = atof(argv[8]);

    double yc = atof(argv[9]);
    double sig_y = atof(argv[10]);
    double py = atof(argv[11]);


    if (instruction == "probability") {

        // Makes filename
        string filename = make_filename(instruction, n_slits);

        // Makes the box
        Box box = Box(n_slits,h, dt, xc, yc, px, py, sig_x, sig_y, barrier_potential);


        // Evolve for the duration and writes total probability to file for each step
        box.evolve_probability(time, filename);

    }
}