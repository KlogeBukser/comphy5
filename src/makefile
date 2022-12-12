
all: all-c clean

all-mac: compile-mac probability clean

all-c: compile probability

heat: compile heat_sim heatmap
interference: compile heat_sim interference_plot
evo: compile evo_sim animate

compile:
	g++ -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo -fopenmp

compile-mac:
	g++-12 -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo -fopenmp

probability: # instruction, n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, V0, T
	main.out probability 0 0.005 0.000025 0.25 0.05 200 0.5 0.05 0 0 0.008
	main.out probability 2 0.005 0.000025 0.25 0.05 200 0.5 0.1 0 1e10 0.008
	python3 time_prob.py


heat_sim: # instruction, n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, V0, T[]
	main.out heat 2 0.005 0.000025 0.25 0.05 200 0.5 0.2 0 1e5 0 0.001 0.002

evo_sim: # instruction, n_slits, h, dt, xc, sig_x, px, yc, sig_y, py, V0, T
	main.out evo 2 0.005 0.000025 0.25 0.05 200 0.5 0.2 0 1e5 0.002

animate:
	python3 animate.py

heatmap:
	python3 heatmap.py

interference_plot:
	python3 interference.py

clean:
	rm f *.out textfiles/*.txt *.txt

clean-plots:
	-rm -f plots/*.pdf *pdf