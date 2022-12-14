
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

infile = open("textfiles/heat_2.txt")
lines = infile.readlines()
infile.close()

n_values = int(lines[0]) + 2

for line in lines[1:]:
    time, data = line.split(':')
    time = 1000*float(time)
    value_sets = data.split(',')

    prob = np.zeros((n_values,n_values))
    real = np.zeros((n_values,n_values))
    imag = np.zeros((n_values,n_values))

    for values in value_sets:
        values = values.split()
        j = int(values[0]) + 1
        i = int(values[1]) + 1

        prob[i,j] = float(values[2])
        real[i,j] = float(values[3])
        imag[i,j] = float(values[4])    

    # Some settings
    fontsize = 12
    x_min, x_max = 0, 1
    y_min, y_max = 0, 1
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(prob))


    for elem,name,title_part in zip([prob,np.sqrt(prob),real,imag],["prob","root","real","imag"], ["P", "Square root of p","Real part of p", "Imaginary part of p"]):
       
        title = title_part + "robability distribution after " + str(time) + " ms"
    
        # Create figure
        fig = plt.figure()
        ax = plt.gca()
        
        img = ax.imshow(elem, extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm = norm)

        # Axis labels
        plt.xlabel("x", fontsize=fontsize)
        plt.ylabel("y", fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)

        # Add a colourbar
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label("Density", fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)

        plt.title(title)

        plt.savefig("plots/distribution_" + name + str(time) + ".pdf")
        plt.close()

