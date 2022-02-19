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
import os
# %%  BULK

os.system("rm -r bands.out")

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

maxCount = 20 * 8
count = 0


def writeBasicStuff(file):
    file.write(
        """#                                                                     
# :   :::   :::     ::::    ::::  :::::::::   ::::::::                
# :+:   :+: :+: :+:   +:+:+: :+:+:+ :+:    :+: :+:    :+              
#  +:+ +:+ +:+   +:+  +:+ +:+:+ +:+ +:+    +:+ +:+    +:+             
#   +#++: +#++:++#++: +#+  +:+  +#+ +#++:++#+  +#+    +:+             
#    +#+  +#+     +#+ +#+       +#+ +#+    +#+ +#+    +#+             
#    #+#  #+#     #+# #+#       #+# #+#    #+# #+#    #+#             
#    ###  ###     ### ###       ### #########   ########              
#                                                                     
#                                                                     
#       Version 5.0.4 Revision 19595 Hash 896bffc02                   
#                        Branch is                                    
#                    MPI+HDF5_IO Build                                
#                http://www.yambo-code.org                            
#
electrons                        # [R] Electronic properties
bnds                             # [R] Bands
INTERP_mode= "BOLTZ"                # Interpolation mode (NN=nearest point, BOLTZ=boltztrap aproach)
INTERP_Shell_Fac= 20.00000       # The bigger it is a higher number of shells is used
cooIn= "rlu"                     # Points coordinates (in) cc/rlu/iku/alat
cooOut= "rlu"                    # Points coordinates (out) cc/rlu/iku/alat
CIRCUIT_E_DB_path= "none"        # SAVE obtained from the QE `bands` run (alternative to %BANDS_kpts)
BANDS_path= ""                   # High-Symmetry points labels (G,M,K,L...) also using composed positions (0.5xY+0.5xL).
GfnQPdb= "E < ./run_all_kpoints/ndb.QP" # read the precalculated QP corrections
BANDS_steps=30
% BANDS_bands
""")


for bands in [
        "80 | 95 |", "96 | 110 |", "111 | 125 |", "126 | 140 |", "141 | 155|",
        "156 | 170 |", "171 | 185 |", "186 | 200 |", "201 | 215 |", "216 | 230 |", "231| 245 |", "246 | 260 |", "261 | 275 |"
]:
    for depthFrac in np.linspace(0, 1, 20):
        with open(f"ypp_bands.in", "w") as file:
            writeBasicStuff(file)
            file.write(bands + "\n")
            file.write("%\n")
            file.write("%BANDS_kpts \n")
            for frac in tenStepsNoEnd:
                kpoint = (1 -
                          frac) * start + middle * (frac) + depthFrac * (D - B)
                file.write(
                    f"{kpoint[0]:12.8} |{kpoint[1]:12.8} |{kpoint[2]:12.8}|\n")
            for frac in tenSteps:
                kpoint = (1 - frac) * middle + end * (frac) + depthFrac * (D -
                                                                           B)
                file.write(
                    f"{kpoint[0]:12.8} |{kpoint[1]:12.8} |{kpoint[2]:12.8}|\n")
            file.write("%\n")
        count += 1
        print(f"Doing iteration {count} of {maxCount}")
        os.system(f"ypp -F ypp_bands.in -C bands.out")
