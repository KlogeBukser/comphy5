
import numpy as np
import matplotlib.pyplot as plt

for n_slits in ['0','2']:
	infile = open("textfiles/probability_" + n_slits + ".txt")
	lines = infile.readlines()
	infile.close()

	n_values = len(lines)
	time = np.empty(n_values)
	probability = np.empty(n_values)

	for i in range(n_values):
		value_set = lines[i].split(',')
		time[i] = 1000*float(value_set[0])
		probability[i] = 1 - float(value_set[1])

	plt.plot(time,probability, label = n_slits + " slits")

plt.title("Deviation of the total probability")
plt.xlabel("Time (ms)")
plt.ylabel("Deviation")
plt.legend()
plt.savefig("plots/tot_prob.png")
plt.show()