A parallel N-Bodies solver, built using fortran-package-manager. To build, first install fortran package manager, and nvhpc. Then: 
```
chmod +x ./build_nv.sh 
./build_nv.sh 
```
Will generate 4 executables: 
1. `mynbodies_serial_out` : A basic serial N-body solver
2. `mynbodies_myomp_out` : An OpenMP parallelized solver for thread based CPU parallelization
3. `mynbodies_omp_out` : An LLM generated OpenMP solver. Similar performance was seen.
4. `mynbodies_acc_out` : An OpenACC solver for GPU parallelization, featuring a one-time data transfer avoiding redundant and repeated host-device copies.

These will generate .csv file outputs which can be visualized using Paraview. I also experimented with LAMMPS dump format but dumped that idea as ParaView can be a bit glitchy. 
Shout out to https://github.com/lele394 for help with the pre-processing scripts. Julia pre-processing scripts in the Input folder were LLM generated.
