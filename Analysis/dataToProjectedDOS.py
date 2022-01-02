#!/usr/bin/env python3

#%%
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LogNorm

#############################################
# INPUT PARAMETERS
#############################################
energy_range = np.array([0, 15])
energyShift = -9.869
outputFileName = "output.pdf"

#############################################
# FUNCTION DEFINITIONS
#############################################


def energiesToDensity(energy_range: np.ndarray,
                      allEnergiesAtKPoint: np.ndarray,
                      resolution: int = 100,
                      sigma: float = 0.1):
    """ Convolution of an array of energies with a gaussian of width sigma. 
    
    In this script this sigma is the spreading (or smearing) along the energy direction of the plot.

    The energy_range consists of 2 numbers that indicate the range over which the convolution is taken.

    The resolution parameter is used to get more data points along the energy range.
    """
    x_range = np.linspace(energy_range[0], energy_range[1], resolution)
    diff = allEnergiesAtKPoint[:, np.newaxis] - x_range
    convolvedValues = 1.0 / (np.sqrt(2 * np.pi * np.power(sigma, 2))) * np.exp(
        -np.power(diff, 2) / (2 * np.power(sigma, 2)))
    res = np.sum(convolvedValues, axis=0)
    res = res / np.linalg.norm(res)
    return res


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


def processDFT(orderedFileList: list,
               energy_range: np.ndarray,
               convFactor: float = 13.6,
               shift: float = 0.0,
               resolution=100,
               sigma=0.1) -> np.ndarray:
    """ processes the res_<nr> output files. """
    results = []
    for file in orderedFileList:
        data = np.loadtxt(file)
        energies = convFactor * data[:, 1] + shift
        yval = energiesToDensity(energy_range, energies, resolution, sigma)
        results.append(yval)
    return np.array(results).transpose()


def visualizeNumExp(ax,
                    basename: str,
                    shift: float,
                    addBar: bool = False,
                    sigma: float = 0.1,
                    Ef: float = None) -> None:
    """ Visualize the data from the numerical experiments as a imshow graphic.

    This function is not general, it assumes a K-path of 59 steps.
    """
    fileList = []
    for i in range(0, 59):
        fileList.append(f"../NumericalExperiments/{basename}/res_{i}.dat")
    results = processDFT(fileList,
                         energy_range,
                         shift=energyShift + shift,
                         resolution=400,
                         sigma=sigma)

    c = ax.imshow(-np.power(results, 1),
                  cmap="viridis",
                  extent=[0, 1, energy_range[0], energy_range[1]],
                  aspect="auto",
                  interpolation="bilinear",
                  origin="lower")
    ax.set_ylim(energy_range)
    ax.set_title(f"{basename}, shift {energyShift+shift}eV")
    ax.set_xticks(tick_loc)
    ax.set_xticklabels(['R', 'G', 'V'])
    if addBar:
        fig.colorbar(c, ax=ax)
    if not Ef is None:
        ax.axhline(Ef + energyShift + shift, color='red', lw='2')

# %%

##############################################
# MAIN 1: Making a plot of a numerical experiment
##############################################

# Get the Experimental Data ...
GtoV = np.loadtxt("../RawData/ARRES1_G-to-V.txt")
RtoG = np.loadtxt("../RawData/ARRES1_R-to-G.txt")
arres1 = np.concatenate([RtoG, GtoV], axis=1)
symmPoints1 = np.array([RtoG.shape[1], GtoV.shape[1]])
# ... and get the tick locations
tick_loc = np.concatenate([[0], np.cumsum(symmPoints1)]) / np.sum(symmPoints1)

# Create reference plot based of arres data
arres1 = arres1 / np.max(arres1)
fig, ax = plt.subplots(1, 2, figsize=(15, 7.5))
visualizeARRES(fig,
               ax[0],
               arres1,
               x_labels=['R', 'G', 'V'],
               label_loc=tick_loc,
               cmap='viridis',
               addBar=True)
ax[0].set_ylim([0, 15])

# Create a plot of a numerical experiment 
# choices are: pent_1layer, pent_2layer, pent_3layer, pent_4layer, pent_bulk

# NB calculations for this plot might take a while, the convolutions are not
# optimized in any way
visualizeNumExp(ax[1], "pent_bulk", 0.0, addBar=True, sigma=0.2)

plt.savefig(outputFileName, dpi=150)

# This will crop the whitespace of the figure, very useful, but it will only
# work on a linux distro with pdfcrop available.
# os.system(f"pdfcrop --margin 5 {outputFileName } {outputFileName }")

# %%
##############################################
# MAIN 2: Making many plots and combining them
##############################################
energy_range = np.array([-20, 15])

fig, ax = plt.subplots(1, 6, figsize=(38, 12))
visualizeARRES(fig,
            ax[0],
            arres1,
            x_labels=['R', 'G', 'V'],
            label_loc=tick_loc,
            cmap='viridis')
ax[0].set_ylim(energy_range)

# Create plots for the numerical experiments
# NB calculations for this plot might take a while, the convolutions are not
# optimized in any way
visualizeNumExp(ax[1], "pent_1layer",  2.924, Ef=0.045)
visualizeNumExp(ax[2], "pent_2layer",  1.972, Ef=1.013)
visualizeNumExp(ax[3], "pent_3layer",  1.564, Ef=1.50)
visualizeNumExp(ax[4], "pent_4layer", 1.224, Ef=1.78)
visualizeNumExp(ax[5], "pent_bulk",  0.0, Ef=2.95, addBar=True)

plotname = f"allLayers.pdf"
plt.savefig(plotname, dpi=150)

# %%
