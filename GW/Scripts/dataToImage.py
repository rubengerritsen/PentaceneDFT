#!/usr/bin/env python3

# %%
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

# I mount my cube1 home folder on remotes/cube1, for this to work you need to
# provide your own path
pathToFiles = "/Users/bbaumeie/work/PentaceneDFT/GW/bands_DFT.out/"
data_DFT = np.loadtxt(pathToFiles + f"o.bands_interpolated")[:, 1:-3]
for i in range(1, 260):
    if i == 99:  # yeah weird, but I think this one is corrupted or something, so I skip it
        continue
    data_DFT = np.concatenate([data_DFT,
        np.loadtxt(pathToFiles + f"o.bands_interpolated_{i:02d}")[:, 1:-3]
    ], axis=1)


pathToFiles = "/Users/bbaumeie/work/PentaceneDFT/GW/bands_GW.out/"
data_GW = np.loadtxt(pathToFiles + f"o.bands_interpolated")[:, 1:-3]
for i in range(1, 260):
    if i == 99:  # yeah weird, but I think this one is corrupted or something, so I skip it
        continue
    data_GW = np.concatenate([data_GW,
        np.loadtxt(pathToFiles + f"o.bands_interpolated_{i:02d}")[:, 1:-3]
    ], axis=1)    

def energiesToDensity(energy_range: np.ndarray,
                      allEnergiesAtKPoint: np.ndarray,
                      resolution: int = 100,
                      sigma: float = 0.03):
    """ Convolves an array of energies with a gaussian of width sigma. """
    x_range = np.linspace(energy_range[0], energy_range[1], resolution)
    diff = allEnergiesAtKPoint[:, np.newaxis] - x_range
    convolvedValues = 1.0 / (np.sqrt(2 * np.pi * np.power(sigma, 2))) * np.exp(
        -np.power(diff, 2) / (2 * np.power(sigma, 2)))
    res = np.sum(convolvedValues, axis=0)
    return res

energy_range = np.array([0, 12])
# after Yambo, the interpolated bands are energy shifted such that 
# EFermi is in the middle of the Gap.
energyShift = -3.5  # manual shift for the DFT calculation = -Work function


fig, ax = plt.subplots(1, 3, figsize=(21, 7.5))
results = []
for i in range(61):
    yval = energiesToDensity(energy_range,
                              data_DFT[i,:] + energyShift,
                              sigma=0.2,
                              resolution=100)
    results.append(yval)

results = np.array(results).transpose()

c = ax[0].imshow(-np.power(results,1),
              cmap="viridis",
              extent=[0, 1, energy_range[0], energy_range[1]],
              aspect="auto",
              interpolation="bilinear",
              origin="lower")
ax[0].set_ylabel("Energy (eV)")
ax[0].set_title("DFT")

# Experimental data, this comes from the RawData folder of the repo

def visualizeARRES(fig,
                   ax,
                   arresData,
                   x_labels=[0, 1],
                   label_loc=[0, 1],
                   extent=[0, 1, -3, 25],
                   cmap="terrain",
                   addBar=False):
    """ Visualizes the ARRES data with an imshow graphic.

    Parameters:
    fig       (figure ) : matplotlib figure
    ax        (subplot) : a subplot of a matplotlib figure.
    arresData (np.array): Matrix with ARRES intensities (k-path along x-axis),
                          assumed to be normalized (values between 0 and 1).
    x_labels  (list str): Labels of K-points (G, R etc.)
    label_loc (list flt): list of x label locations
    extent     (4-array): [min-x, max-x, min-y, max-y].
    """
    c = ax.imshow(
        arresData,
        cmap=cmap,
        norm=LogNorm(vmin=0.01, vmax=1),
        extent=extent,
        aspect="auto",
        interpolation="nearest",
        origin="lower",
    )
    ax.set_ylabel("Energy (eV)")
    ax.set_xticks(label_loc)
    ax.set_xticklabels(x_labels)
    ax.set_title("ARRES DATA")
    if addBar:
        fig.colorbar(c, ax=ax)

GtoV = np.loadtxt("../../RawData/ARRES1_G-to-V.txt")
RtoG = np.loadtxt("../../RawData/ARRES1_R-to-G.txt")
arres1 = np.concatenate([RtoG, GtoV], axis=1)
symmPoints1 = np.array([RtoG.shape[1], GtoV.shape[1]])
# ... get tick locations ...
tick_loc = np.concatenate([[0], np.cumsum(symmPoints1)]) / np.sum(symmPoints1)
# ... and normalize.

# CREATE REFERENCE PLOT, 1 LAYER ARRES
arres1 = arres1 / np.max(arres1)
visualizeARRES(fig,
            ax[1],
            arres1,
            x_labels=['R', 'G', 'V'],
            label_loc=tick_loc,
            cmap='viridis',
            addBar=False)
ax[1].set_ylim([energy_range[0], energy_range[1]])


# And a reference GW calculation
# after Yambo, the interpolated bands are energy shifted such that 
# EFermi is in the middle of the Gap.
energyShift = -3.87  # manual shift for the GW calculation = -Work function + shifted Efermi

results = []
for i in range(61):
    yval = energiesToDensity(energy_range,
                              data_GW[i,:] + energyShift,
                              sigma=0.2,
                              resolution=100)
    results.append(yval)

results = np.array(results).transpose()

c = ax[2].imshow(-np.power(results,1),
              cmap="viridis",
              extent=[0, 1, energy_range[0], energy_range[1]],
              aspect="auto",
              interpolation="bilinear",
              origin="lower")
ax[2].set_ylabel("Energy (eV)")
ax[2].set_title("GW")

# %%
