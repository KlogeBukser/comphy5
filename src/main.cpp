
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <armadillo>

#include "headers/free_functions.h"
#include "headers/Box.h"



constexpr auto ERR_MISSING_PARAMS = 1;

// Standard library
using namespace std;
// Armadillo
using namespace arma;
// Free functions
using namespace free_funcs;

int main(int argc, char* argv[])
{
    if (!eval_argc(argc)) {
        return ERR_MISSING_PARAMS;
    }

    std::string instruction = argv[1];
    int n_slits = atoi(argv[2]);

    double h = atof(argv[3]);
    double dt = atof(argv[4]);

    double xc = atof(argv[5]);
    double sig_x = atof(argv[6]);
    double px = atof(argv[7]);

    double yc = atof(argv[8]);
    double sig_y = atof(argv[9]);
    double py = atof(argv[10]);
    double barrier_potential = atof(argv[11]);

    if (instruction == "probability") {

        double time = atof(argv[12]);

        // Makes filename
        string filename = make_filename(instruction, n_slits);

        // Makes the box
        Box box = Box(n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, barrier_potential);


        // Evolve for the duration and writes total probability to file for each step
        box.write_probability(time, filename);

    }

    if (instruction == "heat") {

        // Makes filename
        string filename = make_filename(instruction, n_slits);

        // Makes the box
        Box box = Box(n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, barrier_potential);

        vec times(argc - 12);
        for (int i = 0; i < argc - 12; i++) {
            times[i] = atof(argv[i+12]);
        }

        // Evolve for the duration and writes total probability to file for each step
        box.write_state(times, filename);

    }

    if (instruction == "evo") {

        double time = atof(argv[12]);

        // Makes filename
        string filename = make_filename(instruction, n_slits);

        // Makes the box
        Box box = Box(n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, barrier_potential);

        // Evolve for the duration and writes total probability to file for each step
        box.write_state(time, filename);

    }
}