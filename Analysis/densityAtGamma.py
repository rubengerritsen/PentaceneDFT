#%%

import numpy as np
import matplotlib.pyplot as plt
import os


def energiesToDensity(energy_range: np.ndarray,
                      allEnergiesAtKPoint: np.ndarray,
                      resolution: int = 100,
                      sigma: float = 0.1):
    """ Convolves an array of energies with a gaussian of width sigma. """
    x_range = np.linspace(energy_range[0], energy_range[1], resolution)
    diff = allEnergiesAtKPoint[:, np.newaxis] - x_range
    convolvedValues = 1.0 / (np.sqrt(2 * np.pi * np.power(sigma, 2))) * np.exp(
        -np.power(diff, 2) / (2 * np.power(sigma, 2)))
    res = np.sum(convolvedValues, axis=0)
    return res


# point 29 is the gamma point
bulk = np.loadtxt("../NumericalExperiments/pent_bulk/res_29.dat")
layer1 = np.loadtxt("../NumericalExperiments/pent_1layer/res_29.dat")
layer2 = np.loadtxt("../NumericalExperiments/pent_2layer/res_29.dat")
layer3 = np.loadtxt("../NumericalExperiments/pent_3layer/res_29.dat")
layer4 = np.loadtxt("../NumericalExperiments/pent_4layer/res_29.dat")

energy_range = np.array([-10, 20])
res = 300
x = np.linspace(energy_range[0], energy_range[1], res)

convFactor = 13.6

global_shift = -9.869

layer1_y = energiesToDensity(energy_range,
                             convFactor * layer1[:, 1] + 2.924 + global_shift,
                             resolution=res,
                             sigma=0.3)
layer2_y = energiesToDensity(energy_range,
                             convFactor * layer2[:, 1] + 1.972 + global_shift,
                             resolution=res,
                             sigma=0.3)
layer3_y = energiesToDensity(energy_range,
                             convFactor * layer3[:, 1] + 1.564 + global_shift,
                             resolution=res,
                             sigma=0.3)
layer4_y = energiesToDensity(energy_range,
                             convFactor * layer3[:, 1] + 1.224 + global_shift,
                             resolution=res,
                             sigma=0.3)
bulk_y = energiesToDensity(energy_range,
                           convFactor * bulk[:, 1] + global_shift,
                           resolution=res,
                           sigma=0.3)
fig, ax = plt.subplots(figsize=(10, 6))

# offset is used to shift densities relative to eachother so they no longer overlap
offset = 0
ax.plot(x, layer1_y / 2 + 0 * offset, ls=':', lw=2, label="1layer")
ax.plot(x, layer2_y / 2 + 1 * offset, ls=':', lw=2, label="2layer")
ax.plot(x, layer3_y / 3 + 2 * offset, ls=':', lw=2, label="3layer")
ax.plot(x, layer4_y / 3.4 + 3 * offset, ls=':', lw=2, label="4layer")
ax.plot(x, bulk_y + 4 * offset, ls='--', lw=2, label="bulk")
ax.legend()
fig.show()

plotname = "densitiesAtGamma.pdf"
plt.savefig(plotname, dpi=150)
os.system(f"pdfcrop --margin 5 {plotname} {plotname}")

# %% Compare the measured cut to some calculation (in this case the 1 layer calculation)
fig, ax = plt.subplots(figsize=(10, 6))


iv_energies = np.load("../RawData/IV_energies.npy")
iv_intensities = np.load("../RawData/IV_intensities.npy")
ax.plot(x, layer1_y / 2, ls=':', lw=2, label="1layer")
ax2 = ax.twinx()
ax2.plot(iv_energies,
         -np.log(iv_intensities) - 0.45*iv_energies,
         ls='--',
         lw=1,
         marker='s',
         markersize=2,
         color='red',
         label="IV curve")
ax2.set_xlim(energy_range)
ax2.legend(loc='upper left')
plotname = f"densityAtGamma_tmp.pdf"
plt.savefig(plotname, dpi=150)
os.system(f"pdfcrop --margin 5 {plotname} {plotname}")

# %%
