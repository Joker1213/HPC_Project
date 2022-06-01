static char help[] = "Solves The Transient Heat Problem By Explicit Method.\n\n";

#include <petscksp.h>
#include <math.h>

#define pi acos(-1)

int main(int argc, char **args){
    Vec u,u_new,b;
    Mat A;
    PetscInt i, j, col[3], rstart, rend, nlocal, its;
    PetscInt n=100;
    PetscReal L, dx, dt;
    PetscReal kappa, rho, c, CFL, value[3];
    PetscScalar u0;
    PetscInt size, rank;
    PetscErrorCode ierr;
    L = 1.0;dx = L/n;dt = 0.00001;its = (int)(1/dt);
    kappa = 1.0;rho = 1.0;c = 1.0;
    CFL = kappa*dt/(rho*c*dx*dx);

    ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "CFL = %f\n", CFL);CHKERRQ(ierr);
 
    // 创建速度向量u
    ierr = VecCreate(PETSC_COMM_WORLD, &u);CHKERRQ(ierr);
    ierr = VecSetSizes(u,PETSC_DECIDE,n+1);CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);CHKERRQ(ierr);
    ierr = VecDuplicate(u, &u_new);CHKERRQ(ierr);
    ierr = VecDuplicate(u, &b);CHKERRQ(ierr);

    // 速度初始化
    for (j = 0; j <= n; j++)
    {
        if(j==0 || j == n){
            u0 = 0.0;
            ierr = VecSetValue(u, j, u0, INSERT_VALUES);CHKERRQ(ierr);
        }
        else{
            u0 = exp(j*dx);     // u(t=0)=e^x
            ierr = VecSetValue(u, j, u0, INSERT_VALUES);CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
    ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
 
    // 确定向量布局
    ierr = VecGetOwnershipRange(u,&rstart,&rend);CHKERRQ(ierr);
    ierr = VecGetLocalSize(u,&nlocal);CHKERRQ(ierr);

    // 创建系数矩阵A
    ierr = MatCreate(PETSC_COMM_WORLD, &A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,nlocal,nlocal,n+1,n+1);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);
 
    // 系数矩阵初始化
    if (!rstart) 
    {
        rstart = 1;
        i      = 0; col[0] = 0; col[1] = 1; value[0] = 0; value[1] = 0;
        ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    if (rend == n+1) 
    {
        rend = n;
        i    = n; col[0] = n-2; col[1] = n-1; value[0] = 0; value[1] = 0;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    value[0] = CFL; value[1] = 1-2.0*CFL; value[2] = CFL;
    for (i = rstart; i < rend; i++)
    {
        col[0] = i-1;col[1] = i;col[2] = i+1;
        ierr = MatSetValues(A, 1, &i, 3, col, value, INSERT_VALUES);CHKERRQ(ierr);
    }
    
    /* Assemble the martix */
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    /* 设置向量b */
    for (j = 0; j <= n; j++)
    {
        if(j==0 || j == n){
            u0 = 0.0;
            ierr = VecSetValue(b, j, u0, INSERT_VALUES);CHKERRQ(ierr);
        }
        else{
            u0 = dt*sin(L*pi*j*dx)/(rho*c);     
            ierr = VecSetValue(b, j, u0, INSERT_VALUES);CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
    // ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    /* u_new = A*u+b */
    for (i = 0; i < its; i++)
    {
        ierr = MatMultAdd(A, u, b, u_new);CHKERRQ(ierr);        // u_new = A*u+b
        ierr = VecCopy(u_new,u);CHKERRQ(ierr);
    }

    ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    
      
    ierr = VecDestroy(&u);CHKERRQ(ierr); 
    ierr = VecDestroy(&u_new);CHKERRQ(ierr); 
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);

    ierr = PetscFinalize();  

    return ierr;
}