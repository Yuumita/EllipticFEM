#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include "elliptic_solver.hpp"
#include "types.h"
#include "utils.hpp"
#include "mesh.hpp"

#include <math.h>

const int D = 4;
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

long double l2_error(Mesh<long double, D> &mesh, VectorX<long double> &solution) {
    long double ret = 0.0;

    #pragma omp parallel for
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

    std::cout << "Give the amount of divisions of the unit cube: ";
    std::cout.flush();
    int divs; std::cin >> divs;

    std::cout << "Creating mesh... ";
    std::cout.flush();
    Mesh<long double, D> mesh(Mesh<long double, D>::get_unit_cube_triangulation(divs));
    std::cout << "finished!\nTotal elements: " << mesh.get_elements_size() << std::endl;

    EllipticSolver<long double, D> es(&mesh, B, L);

    std::cout << "Assembling and solving system... ";
    std::cout.flush();
    VectorX<long double> u = es.assemble_and_solve();
    std::cout << "finished!" << std::endl;

    std::cerr << u << std::endl;

    long double error = l2_error(mesh, u);
    std::cout << "L2 error:                \t" << error << std::endl;

    u.setZero();
    error = l2_error(mesh, u);
    std::cout << "L2 norm of the solution:\t" << error << std::endl;

    return 0;
}

