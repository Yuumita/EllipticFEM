#include <iostream>
#include <eigen3/Eigen/Dense>

#include "elliptic_solver.hpp"
#include "types.h"
#include "utils.hpp"
#include "mesh.hpp"

const int D = 3;

// -d^2u/dx^2 = f_rhs
double f_rhs(Vector<double, D> coords) {
    double ret = 1;
	for (int i = 0; i < D; i++) {
        ret *= coords[i];
	}
	return ret;
}



int main() {
    LinearForm<double, D> L(f_rhs);
    BilinearForm<double, D> B;

    Mesh<double, D> mesh(Mesh<double, D>::get_unit_cube_triangulation());
    EllipticSolver<double, D> es(&mesh, B, L);

    VectorX<double> u = es.solve();
    return 0;
}

