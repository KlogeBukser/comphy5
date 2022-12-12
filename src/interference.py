
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

infile = open("textfiles/interference.txt")
lines = infile.readlines()
infile.close()

n_values = len(lines)

y_values = np.empty(n_values)
density = np.empty(n_values)


for i in range(n_values):
	data_set = lines[i].split()
	y_values[i] = float(data_set[0])
	density[i] = float(data_set[1])

density = density/(np.linalg.norm(density))


plt.title(r"Probability density at $x = 0.8$ after $t = 0.002$ seconds")
plt.plot(y_values,density)
plt.xlabel("y")
plt.ylabel(r"Probability density $|\psi(x,y,t)|^2$")
plt.show()
	
	
