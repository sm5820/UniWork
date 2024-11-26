import netCDF4 as NC
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import matplotlib
import os
from natsort import natsorted
matplotlib.use('Agg')

# ISINGMODEL_V1
dat_v1 = NC.Dataset("ising_fx.nc", "r", format='NETCDF4')

plotv1 = plt.imshow(dat_v1.variables['ISING LATTICE'][:], cmap='Spectral')
cbar = plt.colorbar(plotv1)
cbar.set_label('Spin', rotation=270, labelpad=20)

plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Ising Lattice with fixed boundaries')
plt.savefig('Ising_lattice_fixed_boundaries.png')
plt.legend(['Ising Lattice'], loc='best')
plt.close()

t_v1 = np.arange(1, len(dat_v1.variables['MAGNETISATION'][:])+1)
plt.plot(t_v1, dat_v1.variables['MAGNETISATION'][:], color = 'k')
plt.xlabel('Timesteps')
plt.ylabel('Average Magnetisation')
plt.title('Average Magnetisation of lattice with fixed boundaries')
plt.savefig('magnetisation_fixed_boundaries.png')
plt.close()

#ISINGMODEL_V2
os.chdir('ising_op')
magnetisation_values = []
lattice_data = []

for filename in natsorted(os.listdir('.')):
        dat = NC.Dataset(str(filename), 'r', format='NETCDF4')
        magnetisation_values.append(dat.variables['MAGNETISATION'][:])
        # print(dat.variables['MAGNETISATION'][:])
        lattice_data.append(dat.variables['ISING LATTICE'][:])

lattice_data = np.array(lattice_data)
magnetisation_values = np.array(magnetisation_values[-1])
# print(lattice_data, magnetisation_values)

#Plotting
os.chdir('..')
t_v2 = np.arange(1, len(magnetisation_values)+1)
plt.plot(t_v2-1, magnetisation_values, color = 'blue', linestyle = '--')
plt.xlabel('Timesteps')
plt.ylabel('Average Magnetisation')
plt.title('Average Magnetisation of lattice with periodic boundaries')
plt.savefig('magnetisation_periodic_boundaries.png')
plt.close()

fig, ax = plt.subplots()

def update_frame(i):
    ax.clear() 
    ax.imshow(lattice_data[i], cmap='Spectral', vmin=-1, vmax=1)
    ax.set_title(f"Lattice Evolution - Timestep {i+1}")
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis') 


img = plt.imshow(lattice_data[0], cmap='Spectral')
colorbar = plt.colorbar(img)
colorbar.set_label('Spin', rotation=270, labelpad=20)

ani = animation.FuncAnimation(fig, update_frame, frames=len(lattice_data), interval=500)

ani.save('lattice_evolution_periodic.gif', writer='pillow', fps=2)

plt.close()

# os.rmdir('ising_op')
#Optional command as files in ising_op create problems for next simulation


