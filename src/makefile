
all: all-c clean

all-mac: compile-mac run5 clean

all-c: compile run5


compile:
	g++ -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo -fopenmp

compile-mac:
	g++-12 -std=c++17 -o main.out main.cpp Box.cpp free_functions.cpp -larmadillo -fopenmp

run4:
	main.out test 4 1

run5:
	main.out test 5 1

run6:
	main.out test 6 1


clean:
	rm f *.out textfiles/*.txt *.txt

clean-plots:
	-rm -f plots/*.pdf *pdf