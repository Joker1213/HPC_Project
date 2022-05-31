static char help[] = "Solves The Transient Heat Problem By Explicit Method.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscvec.h>
#include <math.h>


int main(int argc, char **args){
    DM da;
    Vec u,u_new;
    Vec xlocal, glocal;
    PetscInt j, grid=100, L;
    PetscReal kappa, rho, c, CFL, dx, dt, temp;
    PetscInt i, n=10, rstart, rend, nlocal, its, nghost=2, ifrom[2];
    PetscInt size, rank;
    PetscScalar *array, value;
    PetscErrorCode ierr;
    L = 1;dx = L/grid;dt = 0.01;its = 2500;
    kappa = 1.0;rho = 1.0;c = 1.0;
    CFL = kappa*dt/(rho*c*dx*dx);

    ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);

    nlocal = (int)(grid/size);
    if (rank==0){
        ifrom[1] = grid+2*(int)(grid/size)-1;
        ifrom[0] = nlocal;
    }
    else if(rank==size-1){
        ifrom[0] = 0;
        ifrom[1] = rank*nlocal-1;
    }
    else{
        ifrom[0] = (rank+1)*nlocal;
        ifrom[1] = rank*nlocal-1;
    }
    ierr = VecCreateGhost(PETSC_COMM_WORLD, nlocal, PETSC_DECIDE, nghost, ifrom, &u);
    VecDuplicate(u, &u_new);

    ierr = VecGhostGetLocalForm(u, &xlocal);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(u,&rstart,&rend);CHKERRQ(ierr);
    for(i=rstart;i<rend;i++){
        value = exp(dx*(rank*nlocal+i));
        ierr = VecSetValue(u, &i, exp(i*dx), INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(u, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(u, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);

    /*
        Print out each vector, including the ghost padding region.
    */
    ierr = VecGetArray(xlocal, &array);CHKERRQ(ierr);
    for(i=0;i<nlocal+nghost;i++){
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%" PetscInt_FMT " %g\n",i,(double)PetscRealPart(array[i]));
    }
    ierr = VecRestoreArray(xlocal, &array);CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(glocal,&xlocal);CHKERRQ(ierr);

    VecDestroy(u);
    VecDestroy(u_new);
    VecDestroy(xlocal);
    VecDestroy(glocal);

    PetscFinalize();
    
    return 0;
}