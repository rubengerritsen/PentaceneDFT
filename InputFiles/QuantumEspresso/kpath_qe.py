#!/usr/bin/env python3
#%%
# *-------* C
# |       |
# |       |
# |   *B  |
# |       |
# |       |
# *-------* A
# Depth = D

import numpy as np
# %%  BULK

A = np.array([0.49243, 0.47576, 0.70886])
B = np.array([0.0, 0.0, 0.5])
C = np.array([-0.49113, 0.47358, 0.58201])
D = np.array([0.0, 0.0, -0.5])

start = A
middle = B
end = C


tenSteps = np.linspace(0, 1, 2)
tenStepsNoEnd = np.linspace(0, 1, 1, endpoint=False)
counter = 0
with open("k_path_bulk.txt", "w") as file:
  for depthFrac in np.linspace(0,1,20):
    for frac in tenStepsNoEnd :
        kpoint = (1-frac) *start   + middle * (frac)  + depthFrac * (D-B)
        file.write(
            f"{kpoint[0]:12.8} |{kpoint[1]:12.8} |{kpoint[2]:12.8}|\n"
        )
        counter += 1
    for frac in tenSteps :
        kpoint = (1-frac) *middle   + end * (frac) + depthFrac * (D-B)
        file.write(
            f"{kpoint[0]:12.8} |{kpoint[1]:12.8} |{kpoint[2]:12.8}|\n"
        )
        counter += 1
print(f"Count is {counter}")
# %%
