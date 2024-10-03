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

### Example: Solving a Poisson's equation in 4D.

```cpp
#include <iostream>
#include "elliptic_solver.hpp"
#include "types.h"
#include <math.h>

const int D = 4;
const long double PI = 3.1415926535;

// -Δu = f_rhs
long double f_rhs(Vector<long double, D> point) {
    long double ret = (long double)(D) * PI * PI;
    for(int i = 0; i < D; i++) ret *= sin(point[i] * PI);
    return ret;
}

long double analytic_solution(Vector<long double, D> point) {
    long double ret = 1.0;
    for(int i = 0; i < D; i++) ret *= sin(point[i] * PI);
    return ret;
}

long double l2_error(Mesh<long double, D> &mesh, VectorX<long double> &solution) {
    long double ret = 0.0;

    for(Element<long double, D> *e: mesh.get_elements()) {
        for (const auto &qp : e->get_quadrature_points_rand()) {
            const Vector<long double, D>& xq = qp.first;
            long double wq = qp.second;

            long double u_approx = 0.0;
            for (int i = 0; i < e->get_vertices().size(); i++) {
                if(e->get_vertex(i)->is_boundary()) continue;
                u_approx += solution[e->get_vertex(i)->get_index()] * e->get_function(i)(xq);
            }
            long double u_exact = analytic_solution(xq);

            long double diff = u_approx - u_exact;
            ret += wq * diff * diff;
        }
    }

    return ret;
}


int main() {
    LinearForm<long double, D> L(f_rhs);
    BilinearForm<long double, D> B;

    std::cout << "Give the amount of divisions of the unit cube: "; std::cout.flush();
    int divs; std::cin >> divs;

    std::cout << "Creating mesh... ";
    std::cout.flush();
    Mesh<long double, D> mesh(Mesh<long double, D>::get_unit_cube_triangulation(divs));
    std::cout << "finished!\nTotal elements: " << mesh.get_elements_size() << std::endl;

    EllipticSolver<long double, D> es(&mesh, B, L);

    std::cout << "Assembling and solving system... "; std::cout.flush();
    VectorX<long double> u = es.assemble_and_solve();
    std::cout << "finished!" << std::endl;

    std::cerr << u << std::endl;
    long double error = l2_error(mesh, u);
    std::cout << "L2 error:                \t" << error << std::endl;

    return 0;
}
```

## Technical Details

The mathematical background for the FEM approach to solving elliptic PDEs, along with the theoretical details, is well-documented in literature. For an in-depth explanation of the technicalities used in **EllipticFEM** refer to [this article](https://www.algorithmas.org/articles/elliptic_fem.html) which was written in parallel with this implementation. 

In short, we use simplex elements with linear shape functions to approximate the solution. We use a master simplex and compute all other elements, along with their local shape functions, by using the appropriate affine transformation. In computing integrals over elements we mainly use 3 different quadtrature rules: the centroid rule; taking the d+1 midpoints from the vertices to the centroid (and the centroid itself); taking some (k = 35 by default) uniformly random points inside the simplex.
