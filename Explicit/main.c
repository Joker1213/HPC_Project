static char help[] = "Solves The Transient Heat Problem By Explicit Method.\n\n";

#include <petscksp.h>
#include <petscviewerhdf5.h>
#include <math.h>
#include <time.h>

#define pi acos(-1)

int main(int argc, char **args){
    Vec u,u_new,b,temp;
    Mat A;
    PetscInt i, j, col[3], rstart, rend, nlocal, its;
    PetscInt n=5000, index[3], reload=0;
    PetscReal L, dx, dt, u0;
    PetscReal kappa, rho, c, CFL, value[3];
    PetscScalar data[3];
    PetscViewer viewer;
    PetscInt size, rank;
    PetscBool restart = PETSC_FALSE;     // 重启功能标志，初始为FALSE
    PetscErrorCode ierr;

    clock_t start, end;
    double time;
    start = clock();

    ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
    // 从命令行中读取选项参数
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,NULL,NULL);CHKERRQ(ierr); 
    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-restart", &restart, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsEnd();CHKERRQ(ierr);

    // 设置初始参数
    L = 1.0;dx = L/n;dt = 1e-8;its = (int)(1/dt);
    kappa = 1.0;rho = 1.0;c = 1.0;
    CFL = kappa*dt/(rho*c*dx*dx);

    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "CFL = %f\n", CFL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "its = %d\n", its);CHKERRQ(ierr);
 
    // 创建速度向量u
    ierr = VecCreate(PETSC_COMM_WORLD, &u);CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &temp);CHKERRQ(ierr);
    ierr = VecSetSizes(u,PETSC_DECIDE,n+1);CHKERRQ(ierr);
    ierr = VecSetSizes(temp, 3, PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);CHKERRQ(ierr);
    ierr = VecSetFromOptions(temp);CHKERRQ(ierr);
    ierr = VecDuplicate(u, &u_new);CHKERRQ(ierr);
    ierr = VecDuplicate(u, &b);CHKERRQ(ierr);

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
    /* Assemble the vector */
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
    // ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    // 设置是否重启
    if(restart){
        ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "vector.h5", FILE_MODE_READ, &viewer);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)temp, "condition");CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)u, "explicit solution");CHKERRQ(ierr);
        ierr = VecLoad(temp, viewer);CHKERRQ(ierr);
        ierr = VecLoad(u, viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
        index[0] = rank*3;index[1] = rank*3+1;index[2] = rank*3+2;
        ierr = VecGetValues(temp, 3, index, data);CHKERRQ(ierr);
        dx = data[0];dt = data[1];reload = (PetscInt)data[2];
        // 打印每一个CPU的值
        // ierr = PetscPrintf(PETSC_COMM_SELF, "dx = %f, dt = %f, its = %d, rank = %d\n", dx, dt, reload, rank);CHKERRQ(ierr);
    }
    else{
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
        // ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }

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
        i    = n; col[0] = n-1; col[1] = n; value[0] = 0; value[1] = 0;
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
    // ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    // ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "dx=0.01,dt=0.000001.h5", FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    /* u_new = A*u+b */
    for (i = reload; i < its+1; i++)
    {
        ierr = MatMultAdd(A, u, b, u_new);CHKERRQ(ierr);        // u_new = A*u+b
        ierr = VecCopy(u_new,u);CHKERRQ(ierr);
        if((i)%10 == 0){
            index[0] = rank*3;index[1] = rank*3+1;index[2] = rank*3+2;
            data[0] = dx;data[1] = dt;data[2] = i;
            ierr = VecSetValues(temp, 3, index, data,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecAssemblyBegin(temp);CHKERRQ(ierr);    
            ierr = VecAssemblyEnd(temp);CHKERRQ(ierr);   
            ierr = PetscObjectSetName((PetscObject)temp, "condition");CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)u, "explicit solution");CHKERRQ(ierr);
            ierr = VecView(u, viewer);CHKERRQ(ierr);
            ierr = VecView(temp, viewer);CHKERRQ(ierr);
        }
    }

    ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);     // 查看计算最终结果
    // ierr = VecView(temp, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(u, viewer);CHKERRQ(ierr);                        // 将结果输出到HDF5文件
    
    // 计算时间
    end = clock();
    time = (double)(end - start)/CLOCKS_PER_SEC;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "time = %fs\n", time);CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = VecDestroy(&u);CHKERRQ(ierr); 
    ierr = VecDestroy(&u_new);CHKERRQ(ierr); 
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);

    ierr = PetscFinalize(); 

    return ierr;
}