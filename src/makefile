
all: all-c clean

all-mac: compile-mac probability clean

all-c: compile probability


compile:
	g++ -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo -fopenmp

compile-mac:
	g++-12 -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo -fopenmp

probability: # instruction, n_slits, h, dt, T, xc, sig_x, px, yc, sig_y, py, V0
	main.out probability 0 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.05 0 0;
	main.out probability 2 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.1 0 1e10;


clean:
	rm f *.out textfiles/*.txt *.txt

clean-plots:
	-rm -f plots/*.pdf *pdf