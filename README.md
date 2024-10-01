# EllipticFEM

## Overview

**EllipticFEM** is a C++ implementation for solving elliptic partial differential equations (PDEs) using the Finite Element Method (FEM). Namely, it solves the Laplacian boundary value problem which is a generalization of elliptic PDEs. It supports solving such problems in arbitrary dimensions and arbitrary triangulations/meshes rather than being optimized for 2D or 3D meshes. Therefore, in contrast to what is usually implemented, its value lies primarily in theoretical applications and special mathematical constructions.

## Build

To use **EllipticFEM** ensure you have Eigen3 and C++17 installed. To compile the example `main.cpp` file:

```bash
cd EllipticFEM
mkdir build
cd build
cmake ..
make
./EllipticFEM
```

### Example: Solving the Poisson Equation in 2D.

```cpp
const int D = 2;
const long double PI = 3.1415926535;

// -Δu = f_rhs
long double f_rhs(Vector<long double, D> point) {

    long double ret = (long double)(D) * PI * PI;
    for(int i = 0; i < D; i++) {
        ret *= sin(point[i] * PI);
    }
    return ret;
}

long double analytic_solution(Vector<long double, D> point) {

    long double ret = 1.0;
    for(int i = 0; i < D; i++) {
        ret *= sin(point[i] * PI);
    }
    return ret;
}

int main() {
    LinearForm<long double, D> L(f_rhs);
    BilinearForm<long double, D> B;

    int divs = 16;

    Mesh<long double, D> mesh(Mesh<long double, D>::get_2d_unit_cube_triangulation(divs));
    EllipticSolver<long double, D> es(&mesh, B, L);
    VectorX<long double> u = es.solve();

    std::cerr << u << std::endl;

    long double error = l2_error(mesh, u);
    std::cout << "L2 error:  \t" << error << std::endl;

    return 0;
}
```

## Technical Details

The mathematical background for the FEM approach to solving elliptic PDEs, along with the theoretical details, is well-documented in literature. For an in-depth explanation of the technicalities used in **EllipticFEM** refer to [this article](https://www.algorithmas.org/articles/elliptic_fem.html) which was written in parallel with this implementation. 

In short, we use simplex elements with linear shape functions to approximate the solution. We use a master simplex and compute all other elements, along with their local shape functions, by using the appropriate affine transformation. In computing integrals over elements we mainly use 3 different quadtrature rules: the centroid rule; taking the d+1 midpoints from the vertices to the centroid (and the centroid itself); taking some (k = 35 by default) uniformly random points inside the simplex.