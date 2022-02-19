# YAMBO calculations

## Setup 
To run yambo you don't need much. Just compile it and create a module file similar to the following

```text
#%Module
set basedir /home/ruben/Software/yambo
prepend-path PATH $basedir/bin
prepend-path LD_LIBRARY_PATH $basedir/yambo/lib/external
```

The executable is still present in my `/home/ruben/Software` directory on cube1 so you could use that.


## Doing the calculations

### 0. Self-consistent calculations with Quantum Espresso

The first steps are doing two calculations with quantum espresso.

1. A full SCF calculation (`InputFile/pent_scf.in`)
2. A NSCF calculation (`InputFiles/pent_nscf.in`)

To run a calculation using multiple cores (the two 4s need to be scaled down if you use less than 28 threads)

```bash
mpirun -np 28 /home/ruben/Software/qe/bin/pw.x -nt 4 -nd 4 -i <inputFile>.in
```

The calculations have been done and the calculations can be found in `cube1:/data/ruben/yambo_pwscf/PWSCF`.

### 1. Convert QE data to YAMBO data

In whatever folder we have run the NSCF calculation we have a `<projectname>.save` folder. In this case `pent.save`. 

To convert the data

```bash
cd pent.save
p2y            # with yambo loaded into the path
```

You will get a new folder called `SAVE` and that is the folder YAMBO can work with. YAMBO recommends copying this to a clean folder and start doing the yambo stuff from there. 

### 2. Doing the yambo calculations
I moved the yambo calculations to a different location you can find them in `cube1:/home/ruben/pent_yambo`.

Before you can start a calculation on a new `SAVE` directory, you will need to
initialize yambo. Suppose we copied the `SAVE` directory to `myDirectory` than
in `myDirectory` we need to run yambo without any parameters or input, i.e.

```bash
cd myDirectory
yambo
```

This will generate an `r_setup` file with the general setup. In my directories this has already been done.

#### 2.a A GW calculation
To run a GW calculation in parallel you need to set the right environment variables. I have a script for this `GW/Scripts/runYamboParallel.bash`, on cube1 it is often just named `run`. The input file for a basic GW calculation is `GW/Inputfiles/yambo_pent_gw.in`

This calculation has been run and can be found in `cube1:/home/ruben/pent_yambo/YAMBO`

#### 2.b A full GW calculation
YAMBO only does interpolated band calculations, so we do the GW for all points from the QE calculation and then use the Yambo Post Processing (ypp) tool to interpolate the bands.

The full calculation has been done with partially converged parameters (the convergence experiments can be found in `cube1:/home/ruben/pent_yambo/converge` the python script shows what I did, but this was a costly thing, so I was satisfied quickly and decided to just start doing the bands calculation). The calculations are in `cube1:/home/ruben/pent_yambo/bands`. 

In that folder you will find a script `run_parallel` does the GW calculation on all k-points, the input file is also in this repo `GW/InputFiles/yambo_pent_gw_conv_params.in`, but is called `all_kpoints.in` on cube1.

#### 2.c Interpolating the bands
Finally we interpolate the bands with `ypp`, this is a terrible tool that can only interpolate so many bands at a time. To interpolate all the bands I have a python script (very poorly written, but functioning) `GW/Scripts/interpolateBands.py`. On the server it is called `interpolateBands.py` as well. 
Uncomment the line `GfnQPdb= "E < ./run_all_kpoints/ndb.QP"` if you want to interpolate the uncorrected DFT bands.

After running that all the bands are in `bands.out` and we can post process them to make "heat map" images. 

## To keep in mind
My calculations crashed often because I was out of memory. YAMBO's memory usage scales linearly with the number of processes used, so I found it easiest to just reduce that, 14 threads on cube1 worked for me.

## Making the picture

There is a python file `Scripts/dataToImage.py` that shows you how to convert the data to a heat map picture.
