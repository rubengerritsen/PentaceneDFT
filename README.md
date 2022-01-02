# Pentacene Thin Film DFT Data  

## The Data
In the folder `NumericalExperiments` you can find all the data. The different
numerical experiments are named according to the system, hence `pent_1layer` is
a one layer system of pentacene. 

In every experiment folder you will find a list of `.dat` files. They contain
all the energy levels of the different bands. They are ordered as follows, the
number after the underscore in the file name is the location along the K-path
(e.g. R-G-V). Zero means the start and the largest number is the end of the
K-path. Inside the `.dat` file there is also a K-path but only in the z
direction. 

In this way you get a "sheet" of kpoints along the K-path specified and from
this we obtain the projected density of states. Which happens with some analysis
scripts.

## The Analysis
In the folder `Analysis` I have collected some python scripts that help you with
converting the data to pictures. For explanation see the comments in the files.

## Raw Data
In the folder `RawData` is the data from the arres measurements used to make the
reference plot.