#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <vector>
#include "utils.hpp"


#include "vertex.hpp"

#include <Eigen/Dense>
#include "types.h"

/// @brief The element class
/// @tparam Tp Type of coordinates
/// @tparam D The dimension
template <typename Tp, int D>
class Element {
public:
    Element* parent; // points to master element, null if *this is the master element
    std::vector<Vertex<Tp, D>*> vertices;
    std::vector<LinearFunction<Tp>> funcs; // local shape functions
    size_t index;
    AffineTransformation<Tp> transform;
};


template <typename Tp, int D>
Element<Tp, D> get_master_simplex() {
    Element<Tp, D> *m = new Element<Tp, D>;
    m->parent = nullptr;

    m->vertices.resize(D+1);
    m->vertices[0] = new Vertex<Tp, D>(0);
    for(int i = 1; i <= D; i++) {
        Vertex<Tp, D> *vertex = new Vertex<Tp, D>(0);
        (*vertex)[i] = Tp(1);
        m->vertices[i] = vertex;
    }

    std::vector<Tp> c(D+1, Tp(-1));
    c[D] = Tp(1);
    m->funcs.push_back(LinearFunction<Tp>(c));
    c.assign(D+1, 0);
    for(int i = 1; i <= D; i++) {
        c[i] = 1;
        m->funcs.push_back(LinearFunction<Tp>(c));
        c[i] = 0;
    }

    m->transform = AffineTransformation<Tp>(MatrixX<Tp>::Identity(D, D), VectorX<Tp>::Zero(D));
}

template <typename Tp, int D>
Element<Tp, D> get_simplex(const std::vector<Vertex<Tp, D>*> vertices, const Element<Tp, D> &master) {
    Element<Tp, D> *e = new Element<Tp, D>;
    e->parent = &master;

    e->vertices.resize(D+1);
    e->vertices[0] = vertices[0];
    e->transform.b = vertices[0].coords;
    e->transform.A.resize(D, D);
    for(int i = 1; i <= D; i++) {
        e->vertices[i] = vertices[i];
        e->tranform.A.col(i - 1) = vertices[i].coords; 
    }

    AffineTransformation<Tp> transform_inv = e->transform.inverse();
    e->funcs.resize(D+1);
    for(int i = 0; i < D; i++) { 
        e->funcs[i] = compose_affine_linear(transform_inv, master->funcs[i]);
    }

}

#endif /* ELEMENT_HPP */