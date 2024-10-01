#ifndef ELLIPTIC_SOLVER_HPP
#define ELLIPTIC_SOLVER_HPP

#include "mesh.hpp"
#include <eigen3/Eigen/Sparse>

template <typename Tp, int D>
class EllipticSolver {
private:
    BilinearForm<Tp, D> B;
    LinearForm<Tp, D>   L;

    Mesh<Tp, D> *mesh;

    Eigen::SparseMatrix<Tp> stiffness_matrix; 

    VectorX<Tp> rhs;  
    VectorX<Tp> solution;

public:

    EllipticSolver(Mesh<Tp, D>* _mesh) : mesh(_mesh) {}
    EllipticSolver(Mesh<Tp, D>* _mesh, const BilinearForm<Tp, D> &_B, const LinearForm<Tp ,D> &_L) 
        : mesh(_mesh), B(_B), L(_L) {}

    void assemble();
    // void apply_boundary_conditions();
    VectorX<Tp> solve();

    Eigen::SparseMatrix<Tp> get_stiffness() { return stiffness_matrix; };
    VectorX<Tp> get_rhs() { return rhs; };
};


template <typename Tp, int D>
void EllipticSolver<Tp, D>::assemble() {
    int N = mesh->get_inner_vertices_size();
    stiffness_matrix.resize(N, N); // ???
    stiffness_matrix.reserve(N * (D + 1));
    rhs.resize(N); 
    rhs.setZero();

    for(Element<Tp, D>* e : mesh->get_elements()) {
        for(int i = 0; i < e->get_vertices().size(); i++) {
            Vertex<Tp, D> *I = e->get_vertex(i);
            if(I->is_boundary()) continue;
            rhs.coeffRef(I->get_index()) += L(e->get_function(i), *e);
            for(int j = i; j < e->get_vertices().size(); j++) {
                Vertex<Tp, D> *J = e->get_vertex(j);
                if(J->is_boundary()) continue;
#ifdef DEBUG
                std::cout << "Bilinear of " << I->get_index() << " and " << J->get_index() << " = " << B(e->get_function(i), e->get_function(j), *e) << std::endl;
#endif 
                stiffness_matrix.coeffRef(I->get_index(), J->get_index()) += B(e->get_function(i), e->get_function(j), *e);
                if(I != J) 
                    stiffness_matrix.coeffRef(J->get_index(), I->get_index()) += B(e->get_function(j), e->get_function(i), *e);
            }
        }
    }
#ifdef DEBUG
    std::cout << stiffness_matrix << "\n" << rhs << std::endl << std::endl;
#endif 
}

template <typename Tp, int D>
VectorX<Tp> EllipticSolver<Tp, D>::solve() {
    assemble();
    int N = stiffness_matrix.rows();
    solution.resize(N);
    solution.setZero();

    //Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Tp>> lscg;
    //lscg.setTolerance(Tp(1e-6));
    //lscg.compute(stiffness_matrix);

    //if (lscg.info() != Eigen::Success) {
    //    throw std::runtime_error("Decomposition failed");
    //}

    //solution = lscg.solve(rhs);
    //if (lscg.info() != Eigen::Success) {
    //    throw std::runtime_error("Solving failed");
    //}


    Eigen::ConjugateGradient<Eigen::SparseMatrix<Tp>, Eigen::Lower | Eigen::Upper> CGsolver;
    CGsolver.compute(stiffness_matrix);
    solution = CGsolver.solve(rhs);

    // if(CGsolver.info() != Eigen::Success)
    //     throw std::runtime_error("Can't solve linear system");

    // Eigen::SimplicialLLT<Eigen::SparseMatrix<Tp>> LLTsolver;
    // LLTsolver.compute(stiffness_matrix);

    // if(LLTsolver.info() == Eigen::Success) {
    //     solution = LLTsolver.solve(rhs);
    //     if(LLTsolver.info() != Eigen::Success)
    //         throw std::runtime_error("Can't solve linear system");
    // } else {

    // }

    return solution;
}

#endif /* ELLIPTIC_SOLVER_HPP */