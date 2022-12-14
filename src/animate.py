
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Collect simulation data
infile = open("textfiles/evo_2.txt")
lines = infile.readlines()
infile.close()

n_values = int(lines[0]) + 2

# Initial time, and time step length
init_time, first_frame_data = lines[1].split(':')
t_min = float(init_time)

# Densities for the first frame
init_density = np.zeros((n_values,n_values))

for elem in first_frame_data.split(','):
    
    values = elem.split()
    j = int(values[0]) + 1
    i = int(values[1]) + 1

    init_density[i,j] = float(values[2])


fontsize = 12
x_min, x_max = 0, 1
y_min, y_max = 0, 1


# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(init_density))

# Plot the first frame
img = ax.imshow(init_density, extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label(r"$\Psi(x,y,t)$", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)


def animation_density(i):
    # Normalize the colour scale to the current frame?
    
    time, data = lines[i+1].split(':')
    density = np.zeros((n_values,n_values))

    for elem in data.split(','):
    
        values = elem.split()
        j = int(values[0]) + 1
        i = int(values[1]) + 1

        density[i,j] = float(values[2])

    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(density))
    img.set_norm(norm)

    # Update z data
    img.set_data(density)

    # Update the time label
    time_txt.set_text("t = {:.3e}".format(float(time)))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation_density, interval=1, frames=np.arange(0, len(lines[1:]), 2), repeat=False, blit=0)

# Run the animation!
plt.close()

# # Save the animation
anim.save('plots/animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)

