#!/usr/bin/env python3

# %%

import numpy as np
import matplotlib.pyplot as plt

plt.close("all")

bulk = np.loadtxt(
    "../NumericalExperiments/vaccuum_energies/pentacene-GGA_+basis_bulk-potbar.dat")
layer1 = np.loadtxt(
    "../NumericalExperiments/vaccuum_energies/pentacene-GGA_+basis_1Layer-potbar.dat")
layer2 = np.loadtxt(
    "../NumericalExperiments/vaccuum_energies/pentacene-GGA_+basis_2Layer-potbar.dat")
layer3 = np.loadtxt(
    "../NumericalExperiments/vaccuum_energies/pentacene-GGA_+basis_3Layer-potbar.dat")
layer4 = np.loadtxt(
    "../NumericalExperiments/vaccuum_energies/pentacene-GGA_+basis_4Layer-potbar.dat")

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(bulk[:, 0], bulk[:, 1], label="bulk")
ax.plot(layer1[:, 0], layer1[:, 1]+0.215, label="1layer")
ax.plot(layer2[:, 0], layer2[:, 1]+0.145, label="2layer")
ax.plot(layer3[:, 0], layer3[:, 1]+0.115, label="3layer")
ax.plot(layer4[:, 0], layer4[:, 1]+0.090, label="4layer")
ax.legend()

ax2=ax.twinx()
ax2.plot()

# %%
shifts = np.array([0,0.215,0.145,0.115,0.090]) * 13.6
# %%
