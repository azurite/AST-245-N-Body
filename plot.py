import os
import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

basedir_data = "data/"
basedir_plot = ""
metaname = basedir_data + "meta.txt"

files = [
        "data_n2_circular.dat.output",
        "data_n2_elliptic.dat.output",
        "data_n32.dat.output"
    ]

# plot the positions of the integrations
for file in files:
    fullname = basedir_data + file + "_pos.txt"
    metaname = basedir_data + file + "_meta.txt"

    if os.path.exists(fullname):
        dt, T, eps = loadtxt(metaname)
        data = loadtxt(fullname)
        textstr = "\n".join((
            r'Softening = %f' % (eps),
            r'Time Step = %f' % (dt),
            r'Time Interval = [0, %f]' % (T)
        ))

        fig = plt.figure()
        fig.text(0.05, 0.85, textstr, fontsize=12)
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(
            xs=data[0],
            ys=data[1],
            zs=data[2],
            s=0.1,
            label="Planet Trajectory"
        )
        ax.legend()
        fig.savefig(basedir_plot + file + "_pos.png")


# Plot the energies of the integrations
for file in files:
    fullname = basedir_data + file + "_energy.txt"
    metaname = basedir_data + file + "_meta.txt"

    if os.path.exists(fullname):
        dt, T, eps = loadtxt(metaname)
        energy = loadtxt(fullname)
        time = np.arange(0, T, T / len(energy))
        fig, ax = plt.subplots()

        ax.plot(time, energy, "b", label=file)
        ax.set(
            xlabel='time',
            ylabel='|[E(0) - E(t)]/E(0)|',
            title='Relative Energy Error',
        )
        ax.grid()
        ax.legend()
        fig.savefig(basedir_plot + file + "_energy.png")
