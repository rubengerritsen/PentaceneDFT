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
    #res = res / np.linalg.norm(res) # normalize if required
    return res


bulk = np.loadtxt("../numexp/pent_bulk_new/res_29.dat")
layer1 = np.loadtxt("../numexp/pent_1layer/res_29.dat")
layer2 = np.loadtxt("../numexp/pent_2layer/res_29.dat")
layer3 = np.loadtxt("../numexp/pent_3layer/res_29.dat")
layer4 = np.loadtxt("../numexp/pent_4layer/res_29.dat")

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

# offset is used to shift densities relative to eachother for better comparison
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

# %%


iv_energies = np.load("../RawData/IV_energies.npy")
iv_intensities = np.load("../RawData/IV_intensities.npy")
ax2 = ax.twinx()
ax2.plot(iv_energies,
         -(np.log(iv_intensities) + np.exp(0.001 * iv_energies)),
         ls='--',
         lw=1,
         marker='s',
         markersize=2,
         color='red',
         label="IV curve")
ax2.set_xlim(energy_range)
ax2.legend(loc='upper left')
plotname = f"results/densityAtGamma_tmp.pdf"
plt.savefig(plotname, dpi=150)
os.system(f"pdfcrop --margin 5 {plotname} {plotname}")

fig2, ax2 = plt.subplots(figsize=(10, 6))
ax2.scatter(normal, gw)
ax2.set_xlabel("DFT")
ax2.set_ylabel("GW")
ax2.plot(ax2.get_xlim(), ax2.get_ylim(), ls="--", c=".3")
plotname = "GWvsDFT_val.pdf"
plotname = "results/" + plotname
plt.savefig(plotname, dpi=150)
os.system(f"pdfcrop --margin 5 {plotname} {plotname}")

# %%
