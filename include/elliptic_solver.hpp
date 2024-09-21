#ifndef ELLIPTIC_SOLVER_HPP
#define ELLIPTIC_SOLVER_HPP

#include "mesh.hpp"
#include <eigen3/Eigen/Sparse>

template <typename Tp, int D>
class EllipticSolver {
public:
    EllipticSolver(const Mesh<Tp, D>& _mesh) : mesh(_mesh) {}

    void assemble();
    // void apply_boundary_conditions();
    VectorX<Tp> solve();

private:
    BilinearForm<Tp, D> B;
    LinearForm<Tp, D>   L;

    Mesh<Tp, D> *mesh;
    Eigen::SparseMatrix<Tp> stiffness_matrix;
    VectorX<Tp> rhs;  
    VectorX<Tp> solution;

};


template <typename Tp, int D>
void EllipticSolver<Tp, D>::assemble() {
    int N = mesh->get_vertices().size();
    /// stiffness_matrix.reserve( ... );
    rhs.resize(N); 
    rhs.setZero();
    solution.resize(N);
    solution.setZero();
    
    for(const Element<Tp, D>* e : mesh->elements) {
        for(int i = 0; i < e->vertices.size(); i++) {
            Vertex<Tp, D> *I = e->get_vertex(i);
            for(int j = i; j < e->get_vertices().size(); j++) {
                Vertex<Tp, D> *J = e->get_vertex(j);
                stiffness_matrix(I->get_index(), J->get_index()) += B(e->get_function(i), e->get_function(j), e);
                if(I != J) 
                    stiffness_matrix(J->get_index(), I->get_index()) += B(e->get_function(j), e->get_function(i), e);
            }
            rhs(I) += L(e->get_function(i), e);
        }
    }
}

template <typename Tp, int D>
VectorX<Tp> EllipticSolver<Tp, D>::solve() {
    assemble();
    int N = stiffness_matrix.size();
    solution.resize(N);

    Eigen::SimplicialLLT<Eigen::SparseMatrix<Tp>> LLTsolver;
    LLTsolver.compute(stiffness_matrix);

    if(LLTsolver.info() == Eigen::Success) {
        solution = LLTsolver.solve(rhs);
        if(LLTsolver.info() != Eigen::Success)
            throw std::runtime_error("Can't solve linear system");
    } else {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<Tp>, Eigen::Lower | Eigen::Upper> CGsolver;
        CGsolver.compute(stiffness_matrix);
        solution = CGsolver.solve(rhs);

        if(CGsolver.info() != Eigen::Success)
            throw std::runtime_error("Can't solve linear system");
    }

    return solution;
}

#endif /* ELLIPTIC_SOLVER_HPP */