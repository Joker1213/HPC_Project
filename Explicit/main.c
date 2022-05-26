static char help[] = "Solves The Transient Heat Problem By Explicit Method.\n\n";

#include <petscts.h>
#include <petscksp.h>
#include <petscvec.h>
#include <math.h>


int main(int argc, char **args){
    Vec u,u_new;
    PetscInt j, grid, L;
    PetscReal kappa, rho, c, CFL, dx, dt, temp;
    PetscInt i, n=10, rstart, rend, nlocal, rank, its;
    PetscErrorCode ierr;
    L = 1;grid = 100;dx = L/grid;dt = 0.01;its = 2500;
    kappa = 1.0;rho = 1.0;c = 1.0;
    CFL = kappa*dt/(rho*c*dx*dx);

    ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
    ierr = VecSetSizes(u,PETSC_DECIDE,grid);CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);CHKERRQ(ierr);
    ierr = VecDuplicate(u,&u_new);CHKERRQ(ierr);

    ierr = VecGetOwnershipRange(u, &rstart, &rend);CHKERRQ(ierr);
    ierr = VecGetLocalSize(u,&nlocal);CHKERRQ(ierr);

    for (i = 0; i < grid; i++)
    {
        if(i==0 || i==grid-1){
            ierr = VecSetValue(u, &i, 0, INSERT_VALUES);CHKERRQ(ierr);
        }
        ierr = VecSetValue(u, &i, exp(i*dx), INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
    
    for (i = 1; i < its-1; i++)
    {
        temp = 
    }
    
    
    return 0;
}