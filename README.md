# HPC_Project

This the HPC's Course Project

## Usage

1. Clone the repo
2. Build

   ```
   cd Explicit/Implicit
   make
   ```
3. Run
   if you run the code on the TaiYi Server, just use:

   ```
   bsub < ty_script
   ```

   if you run the code on your computer,follow the command, `-n` means the grid number, `-dt` means the time step, `-restart` means the flag of the restart.

   ```
   mpirun -np 5 ./main.out -n 100 -dt=0.00001 -restart=PETSC_FALSE
   ```
