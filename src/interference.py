
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys

n_slits = sys.argv[1]

infile = open("textfiles/interference" + n_slits + ".txt")
lines = infile.readlines()
infile.close()

n_values = len(lines) - 1
x,time = lines[0].split()
time = 1000*float(time)

y_values = np.empty(n_values)
density = np.empty(n_values)


for i in range(n_values):
	data_set = lines[i+1].split()
	y_values[i] = float(data_set[0])
	density[i] = float(data_set[1])

density = density/(np.linalg.norm(density))


plt.title(r"Probability density with " + n_slits + " slit(s) ($x = " + x[:4] + "$, $t = " + str(time) + "$ ms)")
plt.plot(y_values,density)
plt.xlabel("y")
plt.ylabel(r"Probability density $|\psi(x,y,t)|^2$")

plt.savefig("plots/interference" + n_slits + ".png")
plt.close()
	
	
