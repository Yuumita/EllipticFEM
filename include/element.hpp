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
private:
    Element* parent; // points to master element, null if *this is the master element
    std::vector<Vertex<Tp, D>*> vertices;
    std::vector<LinearFunction<Tp, D>> funcs; // local shape functions
    int index;
    AffineTransformation<Tp, D> transform; // affine transformation used to go from the master element to *this

    mutable Vector<Tp, D> centroid;
    mutable Tp volume;
    mutable bool built_data = false;
    static Element* master_simplex;

public:
    Element() {}
    Element(std::vector<Vertex<Tp, D>*> vertices_): vertices(vertices_) {}


    /// @return Volume of the d-simplex given by vertices[0...d] in d-dimensional space.
    Tp get_volume() const {
        if(!built_data) build_data();
        return volume;
    }

    Vector<Tp, D> get_centroid() const {
        if(!built_data) build_data();
        return centroid;
    }

    std::vector<std::pair<Vector<Tp, D>, Tp>> get_quadrature_points() const {
        return get_quadrature_points_2();
    }

    std::vector<std::pair<Vector<Tp, D>, Tp>> get_quadrature_points_1() const {
        return {{get_centroid(), get_volume()}};
    }

    std::vector<std::pair<Vector<Tp, D>, Tp>> get_quadrature_points_2() const {
        std::vector<std::pair<Vector<Tp, D>, Tp>> ret(D + 2);
        Vector<Tp, D> c = get_centroid();
        for(int i = 0; i <= D; i++) {
            Vector<Tp, D> x = (vertices[i]->get_coords() + c) / Tp(2);
            ret[i] = {x, get_volume() / (D + 2)};
        }
        ret[D + 1] = {c, get_volume() / (D + 2)};
        return ret;
    }


    std::vector<std::pair<Vector<Tp, D>, Tp>> get_quadrature_points_rand(int k = 35) const {
        std::vector<std::pair<Vector<Tp, D>, Tp>> ret(k);
        for(int i = 0; i < k; i++) {
            std::vector<Tp> bary(D+1);
            for(int j = 0; j < D; j++) {
                bary[j] = static_cast<Tp>(rand()) / (Tp)(RAND_MAX);
            }
            bary[D] = Tp(1);

            std::sort(bary.begin(), bary.end());

            Vector<Tp, D> point = Vector<Tp, D>::Zero();
            for(int j = 0; j <= D; j++) {
                Tp w = bary[j] - (j ? bary[j - 1] : Tp(0));
                point += w * vertices[j]->get_coords();
            }

            ret[i] = {point, get_volume() / Tp(k)};
        }
        return ret;
    }

    int get_index() { return index; }
    Vertex<Tp, D> *get_vertex(size_t i) { return vertices[i]; }
    LinearFunction<Tp, D> &get_function(size_t i) { return funcs[i]; }

    std::vector<Vertex<Tp, D>*>& get_vertices() { return vertices; }
    std::vector<LinearFunction<Tp, D>>& get_functions() { return funcs; }

    void build_data() const {
        Matrix<Tp, D, D> V(D, D);
        Vector<Tp, D> v0 = vertices[0]->get_coords();
        for(int i = 0; i < D; i++) {
            V.col(i) = vertices[i + 1]->get_coords() - v0;
        }
        volume = std::abs(V.determinant() / factorial<Tp>(D));

        centroid = vertices[0]->get_coords();
        for(int i = 1; i < vertices.size(); i++) {
            centroid += vertices[i]->get_coords();
        }
        centroid /= Tp(D + 1);
        built_data = true;
    }

    void build_data(std::vector<Vertex<Tp, D>*> nvertices) {
        this->vertices = nvertices;
        build_data();
    }


    static Element& get_master_simplex() {
        if(master_simplex != nullptr) 
            return *master_simplex;

        master_simplex = new Element<Tp, D>;
        master_simplex->parent = nullptr;

        master_simplex->vertices.resize(D+1);
        master_simplex->vertices[0] = new Vertex<Tp, D>(Tp(0));
        for(int i = 0; i < D; i++) {
            Vertex<Tp, D> *vertex = new Vertex<Tp, D>(Tp(0));
            vertex->coefRef(i) = Tp(1);
            master_simplex->vertices[i + 1] = vertex;
        }

        Vector<Tp, D> c = Vector<Tp, D>::Constant(D, Tp(-1));
        master_simplex->funcs.push_back(LinearFunction<Tp, D>(c, Tp(1)));
        c = Vector<Tp, D>::Constant(D, Tp(0));
        for(int i = 0; i < D; i++) {
            c[i] = Tp(1);
            master_simplex->funcs.push_back(LinearFunction<Tp, D>(c, Tp(0)));
            c[i] = Tp(0);
        }

        master_simplex->transform = AffineTransformation<Tp, D>(Matrix<Tp, D, D>::Identity(), Vector<Tp, D>::Zero());
        return *master_simplex;
    }


    static Element<Tp, D> get_simplex(const std::vector<Vertex<Tp, D>*> vertices, Element<Tp, D> &master) {
        Element<Tp, D> e(vertices);
        e.parent = &master;

        e.transform.b = e.vertices[0]->get_coords();
        e.transform.A.resize(D, D);
        for(int i = 1; i <= D; i++) {
            e.transform.A.col(i - 1) = e.vertices[i]->get_coords() - e.vertices[0]->get_coords();
        }

#ifdef DEBUG
        std::cout << "Transformation of element between ";
        for(Vertex<Tp, D> *v: e.get_vertices())
            std::cout << *v <<" ";
        std::cout << "\tMatrix A: " << e.transform.A << std::endl;
        std::cout << "\tb: " << e.transform.b << std::endl;
#endif
        if (std::abs(e.transform.A.determinant()) < 1e-8) {
            std::cerr << "Warning: Singular matrix, cannot invert." << std::endl;
            throw std::runtime_error("Singular matrix encountered.");
        }


        AffineTransformation<Tp, D> transform_inv = e.transform.inverse();
        e.funcs.resize(D+1);
        for(int i = 0; i < D; i++) { 
            e.funcs[i] = compose_affine_linear<Tp, D>(transform_inv, master.funcs[i]);
        }


        return e;
    }

    
};


template<typename Tp, int D>
Element<Tp, D>* Element<Tp, D>::master_simplex = nullptr;

#endif /* ELEMENT_HPP */