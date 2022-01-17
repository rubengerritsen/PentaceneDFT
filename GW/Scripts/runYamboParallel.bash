#!/bin/bash

module load yambo

nthreads=1
ncpu=28

export OMP_NUM_THREADS=$nthreads

label=GW_MPI_pentacene
jdir=run_${label}
cdir=run_${label}.out

filein=yambo_gw.in

DIP_CPU="1 $ncpu 1"       # [PARALLEL] CPUs for each role
DIP_ROLEs="k c v"         # [PARALLEL] CPUs roles (k,c,v)
X_CPU="1 1 1 $ncpu 1"     # [PARALLEL] CPUs for each role
X_ROLEs="q g k c v"       # [PARALLEL] CPUs roles (q,g,k,c,v)
X_nCPU_LinAlg_INV=$ncpu   # [PARALLEL] CPUs for Linear Algebra
SE_CPU=" 1 $ncpu 1"       # [PARALLEL] CPUs for each role
SE_ROLEs="q qp b"         # [PARALLEL] CPUs roles (q,qp,b)
SE_Threads=0  

echo "Running on $ncpu MPI threads"
mpirun -np $ncpu yambo -F $filein -J $jdir -C $cdir
