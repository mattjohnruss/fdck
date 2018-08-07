import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from mpl_toolkits.mplot3d import Axes3D
import os
from glob import glob
from sys import argv

file_list = np.sort(glob("./*.csv"))
n_files = len(file_list)
n_points = int(argv[1])
n_var = int(argv[2])
n_timesteps = int(argv[3])
dt = float(argv[4])

data = np.zeros([n_files, n_points, n_var])

for (i, file_path) in enumerate(file_list):
    data[i] = np.genfromtxt(file_path, skip_header=1, usecols=(2, 3, 4, 5, 6))


x = np.linspace(0, 1, n_points)
t = np.arange(0, n_timesteps*dt, dt)

aspect_ratio = x.max()/t.max()
# aspect_ratio = t.max()/x.max()

# X, T = np.meshgrid(x, t)

mpl.rc('mathtext', fontset='cm')

fig = plt.figure()

# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, T, data[0:n_timesteps, :, 1], linewidth=0)

ax = fig.add_subplot(111)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$t$')

im = ax.imshow(data[0:n_timesteps, :, 4], interpolation='none',
               origin='lower', extent=[x.min(), x.max(), t.min(), t.max()],
               aspect=aspect_ratio)

# ax.set_xlabel(r'$t$')
# ax.set_ylabel(r'$x$')

# im = ax.imshow(data[0:n_timesteps, :, 4].T, interpolation='none',
# origin='upper', extent=[t.min(), t.max(), x.min(), x.max()],
# aspect=aspect_ratio)

# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.01)
# plt.colorbar(im, cax=cax)
plt.colorbar(im)

ax.set_title(os.path.basename(os.getcwd()))

# plt.show()
plt.savefig('phi_m.pdf', bbox_inches='tight')
