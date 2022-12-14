
all: all-c clean

all-mac: compile-mac probability heat interference animate clean

all-c: compile execute

execute: probability heat interference animate

compile:
	g++ -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo

compile-mac:
	g++-12 -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo

probability: # instruction, n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, V0, T
	main.out probability 0 0.005 0.000025 0.25 0.05 200 0.5 0.05 0 0 0.008
	main.out probability 2 0.005 0.000025 0.25 0.05 200 0.5 0.1 0 1e10 0.008
	python3 time_prob.py


heat: # instruction, n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, V0, T[]
	main.out heat 2 0.005 0.000025 0.25 0.05 200 0.5 0.2 0 1e5 0 0.001 0.002
	python3 heatmap.py

animate: # instruction, n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, V0, T
	main.out evo 2 0.005 0.000025 0.25 0.05 200 0.5 0.2 0 1e5 0.002
	python3 animate.py


interference:
	main.out heat 1 0.005 0.000025 0.25 0.05 200 0.5 0.2 0 1e5 0.002
	main.out heat 2 0.005 0.000025 0.25 0.05 200 0.5 0.2 0 1e5 0.002
	main.out heat 3 0.005 0.000025 0.25 0.05 200 0.5 0.2 0 1e5 0.002
	python3 interference.py 1
	python3 interference.py 2
	python3 interference.py 3

clean:
	rm f *.out textfiles/*.txt *.txt

clean-plots:
	-rm -f plots/*.pdf *pdf