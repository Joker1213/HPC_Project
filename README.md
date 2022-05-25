# HPC_Project

HPC课程的Project， 关于一维传热方程的数值解(显式方法和隐式方法)

一维传热方程如下：

$$
\begin{aligned}
\rho c \frac{\partial u}{\partial t} - \kappa \frac{\partial^2 u}{\partial x^2} &= \sin(l \pi x)\\
\kappa \frac{\partial u}{\partial x} n_x &= h\\
\kappa &= 1.0\\
u(0,t) &= u(1,t) = 0\\
u|_{t=0} &= e^x
\end{aligned}
$$

边界条件为：


传热方程改写成：

$$
\begin{aligned}
\frac{\partial u}{\partial t} - \frac{\kappa}{\rho c}\frac{\partial^2 u}{\partial x^2} &= \frac{1}{\rho c}\sin(l \pi x)\\
\end{aligned}
$$

## 显式格式(Adams-Bashforth)：

差分格式为：

$$
\left\{ 
    \begin{array}{c}
        \begin{aligned}
        \frac{\partial u}{\partial t} &= \frac{u^{n+1}_j-u^n_j}{\Delta t}\\
        \frac{\partial^2 u}{\partial x^2} &= \frac{u^n_{j+1}-2u^n_j+u^n_{j-1}}{\Delta x^2}
        \end{aligned}
    \end{array}
\right.
$$

$$
\begin{aligned}
\therefore \frac{u^{n+1}_j-u^n_j}{\Delta t} &= \frac{\kappa}{\rho c} \frac{u^n_{j+1}-2u^n_j+u^n_{j-1}}{\Delta x^2} + \frac{1}{\rho c}\sin(l \pi j \Delta x) \\
\implies u^{n+1}_j &= u^n_j + \frac{\kappa \Delta t}{\rho c \Delta x^2} (u^n_{j+1}-2u^n_j+u^n_{j-1})+ \frac{\Delta t}{\rho c}\sin(l \pi j \Delta x)\\
\implies u^{n+1}_j &= \frac{\kappa \Delta t}{\rho c \Delta x^2} u^n_{j+1} + (1-\frac{2 \kappa \Delta t}{\rho c \Delta x^2}) u^n_j + \frac{\kappa \Delta t}{\rho c \Delta x^2}u^n_{j-1} + \frac{\Delta t}{\rho c}\sin(l \pi j \Delta x)\\
make  \frac{\kappa}{\rho c} &= \alpha, \frac{\alpha \Delta t}{\Delta x^2} = \beta = \frac{\kappa \Delta t}{\rho c \Delta x^2}\\
\implies u^{n+1}_j &= \beta u^n_{j+1} + (1-2\beta) u^n_j + \beta u^n_{j-1} + \frac{\Delta t}{\rho c}\sin(l \pi j \Delta x)\\
\end{aligned}
$$

## 隐式格式(Crank-Nicolson)：

差分格式为：

$$
\left\{ 
    \begin{array}{c}
        \begin{aligned}
        \frac{\partial u}{\partial t} &= \frac{u^{n+1}_j-u^n_j}{\Delta t}\\
        \frac{\partial^2 u}{\partial x^2} &= \frac{u^{n+1}_{j+1}-2u^{n+1}_j+u^{n+1}_{j-1}}{\Delta x^2}
        \end{aligned}
    \end{array}
\right.
$$

$$
\begin{aligned}
\therefore \frac{u^{n+1}_j-u^n_j}{\Delta t} &= \frac{\kappa}{\rho c} \frac{u^{n+1}_{j+1}-2u^{n+1}_j+u^{n+1}_{j-1}}{\Delta x^2} + \frac{1}{\rho c}\sin(l \pi j \Delta x)\\
make  \frac{\kappa}{\rho c} &= \alpha, \frac{\alpha \Delta t}{\Delta x^2} = \beta = \frac{\kappa \Delta t}{\rho c \Delta x^2}\\
\end{aligned}
$$

$$
\begin{aligned}
\implies -\beta u^{n+1}_{j+1} + (1+2\beta)u^{n+1}_j - \beta u^{n+1}_{j+1} = u^n_j + \alpha \Delta t \sin(l \pi j \Delta x)
\end{aligned}
$$

解的格式为：

$$
\begin{pmatrix}
    1+2 \beta & -\beta & 0 & 0 & \cdots & 0 & 0 & 0\\
    -\beta & 1+2 \beta & -\beta & 0 & \cdots & 0 & 0 & 0\\
    0 & -\beta & 1+2 \beta & -\beta & \cdots & 0 & 0 & 0\\
    \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots\\
    0 & 0 & 0 & 0 & \cdots & 0 & -\beta & 1+2\beta \\
\end{pmatrix}*
\begin{pmatrix}
    u^{n+1}_2  \\
    u^{n+1}_3 \\
    \vdots\\
    u^{n+1}_{n-1}\\
    u^{n+1}_n\\
\end{pmatrix}=
\begin{pmatrix}
    u^{n}_2 + \alpha \Delta t \sin(l \pi \Delta x) \\
    u^{n}_3 + \alpha \Delta t \sin(l \pi 2 \Delta x)\\
    \vdots\\
    u^{n}_{n-1} + \alpha \Delta t \sin(l \pi (n-2) \Delta x)\\
    u^{n+1}_n + \alpha \Delta t \sin(l \pi (n-1) \Delta x)\\
\end{pmatrix}
$$
