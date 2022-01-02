#!/usr/bin/env python3


import numpy as np
import os

def runBjoernDFT(kpoint: np.ndarray, output: str):
  """ Runs the dft calculation for a given kpoint along the z-axis. """
  zShift = np.array([0.0, 0.0, -0.5])
  kpPlusZ =kpoint+zShift 
  with open("kway_z_scan", "w") as file:
    file.write("begin kway\n")
    file.write("scaling = lattice\n")
    # we use the kpath in the dft calculation to follow the z-axis
    file.write(f"{kpPlusZ[0]:f}d0 {kpPlusZ[1]:f}d0 {kpPlusZ[2]:f}d0\n")
    file.write(f"{kpoint[0]:f}d0 {kpoint[1]:f}d0 {kpoint[2]:f}d0\n")  
    file.write("n_step = 50\n")
    file.write("n_ev_factor = 4\n")
    file.write("end kway\n")
  os.system("make pentacene-bulk-band.dat")
  os.system(f"mv pentacene-bulk-band.dat {output}")


# main points of interest at the top of the brillouin zone
A = np.array([0.49650, 0.00110, 0.62807])
B = np.array([0.0, 0.0, 0.5])
C = np.array([0.00068, 0.49293, 0.80207])

start = A
middle = B
end = C

tenSteps = np.linspace(0, 1, 30)
tenStepsNoEnd = np.linspace(0, 1, 29, endpoint=False)
iter = 0
# trace the line between the first two points of interest
for frac in tenStepsNoEnd:
  kpoint = (1-frac) *start   + middle * (frac)  
  runBjoernDFT(kpoint, f"res_{iter}.dat")
  iter = iter +1
# trace the line between the second and third point
for frac in tenSteps:
  kpoint = (1-frac) *middle   + end * (frac)  
  runBjoernDFT(kpoint, f"res_{iter}.dat")
  iter = iter +1
