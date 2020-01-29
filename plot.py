import os
import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

basedir_data = "data/"
basedir_plot = ""
metaname = basedir_data + "meta.txt"

plotSize = (12,8)

files = [
        "data_n2_circular.dat.output",
        "data_n2_elliptic.dat.output",
        "data_n32.dat.output",
        "data_n64.dat.output",
        "data_n128.dat.output",
        "data_n256.dat.output"
    ]

# plot the positions of the integrations
for file in files:
    fullname = basedir_data + file + "_pos.txt"
    metaname = basedir_data + file + "_meta.txt"

    if os.path.exists(fullname):
        dt, T, eps = loadtxt(metaname)
        data = loadtxt(fullname)
        textstr = "\n".join((
            r'Softening = %g' % (eps),
            r'Time Step = %g' % (dt),
            r'Time Interval = [0, %g]' % (T)
        ))

        fig = plt.figure(figsize=plotSize)
        fig.text(0.05, 0.85, textstr, fontsize=12)
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(
            xs=data[0],
            ys=data[1],
            zs=data[2],
            s=0.02,
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
        fig, ax = plt.subplots(figsize=plotSize)

        labelstr = " ".join((r'dt = %g' % dt, r'${\epsilon}$ = %g' % eps))

        ax.plot(time, energy, "b", label=labelstr)
        ax.set(
            xlabel='time',
            ylabel='|[E(0) - E(t)]/E(0)|',
            title='Relative Energy Error',
        )
        ax.grid()
        ax.legend()
        fig.savefig(basedir_plot + file + "_energy.png")
