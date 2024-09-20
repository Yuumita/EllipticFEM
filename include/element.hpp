#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <vector>
#include "utils.hpp"


#include "vertex.hpp"

#include <eigen3/Eigen/Dense>
#include "types.h"

/// @brief The element class
/// @tparam Tp Type of coordinates
/// @tparam D The dimension
template <typename Tp, int D>
class Element {
public:
    Element() {}
    Element* parent; // points to master element, null if *this is the master element
    std::vector<Vertex<Tp, D>*> vertices;
    std::vector<LinearFunction<Tp>> funcs; // local shape functions
    size_t index;
    AffineTransformation<Tp> transform;

    VectorX<Tp> centroid;
    Tp volume;

    bool built = false;

    /// @return Volume of the d-simplex given by vertices[0...d] in d-dimensional space.
    Tp get_volume() const {
        if(!built) build();
        return volume;
    }

    VectorX<Tp> get_centroid() const {
        if(!built) build();
        return centroid;
    }


    void build() {
        MatrixX<Tp> V(D, D);
        VectorX<Tp> v0 = vertices[0]->get_coords();
        for(int i = 0; i < D; i++) {
            V.col(i) = vertices[i + 1]->get_coords() - v0;
        }
        volume = std::abs(V.determinant() / factorial<Tp>(D));

        centroid = vertices[0]->get_coords();
        for(int i = 1; i < vertices.size(); i++) {
            centroid += vertices[i]->get_coords();
        }
    }

    void build(std::vector<Vertex<Tp, D>*> nvertices) {
        this->vertices = nvertices;
        build();
    }


    static Element* master_simplex = nullptr;
    static Element get_master_simplex() {
        if(master_simplex != nullptr) return *master_simplex;
        Element<Tp, D> *master_simplex = new Element<Tp, D>;
        master_simplex->parent = nullptr;

        master_simplex->vertices.resize(D+1);
        master_simplex->vertices[0] = new Vertex<Tp, D>(0);
        for(int i = 1; i <= D; i++) {
            Vertex<Tp, D> *vertex = new Vertex<Tp, D>(0);
            (*vertex)[i] = Tp(1);
            master_simplex->vertices[i] = vertex;
        }

        std::vector<Tp> c(D+1, Tp(-1));
        c[D] = Tp(1);
        master_simplex->funcs.push_back(LinearFunction<Tp>(c));
        c.assign(D+1, 0);
        for(int i = 1; i <= D; i++) {
            c[i] = 1;
            master_simplex->funcs.push_back(LinearFunction<Tp>(c));
            c[i] = 0;
        }

        master_simplex->transform = AffineTransformation<Tp>(MatrixX<Tp>::Identity(D, D), VectorX<Tp>::Zero(D));
        return *master_simplex;
    }
    
};




template <typename Tp, int D>
Element<Tp, D> get_simplex(const std::vector<Vertex<Tp, D>*> vertices, const Element<Tp, D> &master) {
    Element<Tp, D> e = new Element<Tp, D>;
    e.parent = &master;

    e.vertices = vertices;
    e.transform.b = vertices[0]->get_coords();
    e.transform.A.resize(D, D);
    for(int i = 1; i <= D; i++) {
        e.tranform.A.col(i - 1) = vertices[i]->get_coords();
    }

    AffineTransformation<Tp> transform_inv = e->transform.inverse();
    e.funcs.resize(D+1);
    for(int i = 0; i < D; i++) { 
        e.funcs[i] = compose_affine_linear(transform_inv, master->funcs[i]);
    }

    return e;
}

#endif /* ELEMENT_HPP */