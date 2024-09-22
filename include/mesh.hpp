#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <eigen3/Eigen/Dense>
#include "types.h"
#include "vertex.hpp"
#include "element.hpp"


template <typename Tp, int D>
class Mesh {
private:
    std::vector<Vertex<Tp, D>*>  vertices;
    std::vector<Element<Tp, D>*> elements;

public:
    static Mesh get_unit_cube_triangulation();

    Vertex<Tp, D> *get_vertex(size_t i)   { return vertices[i]; }
    std::vector<Vertex<Tp, D>*> &get_vertices()   { return vertices; }
    size_t get_vertices_size()   { return vertices.size(); }

    Element<Tp, D> *get_element(size_t i) { return elements[i]; }
    std::vector<Element<Tp, D>*> &get_elements()   { return elements; }
    size_t get_elements_size()   { return elements.size(); }

};

template <typename Tp, int D>
Mesh<Tp, D> Mesh<Tp, D>::get_unit_cube_triangulation() {
    Mesh<Tp, D> mesh;

    Vector<Tp, D> coords(D);
    for (int i = 0; i < (1 << D); ++i) {
        for (int j = 0; j < D; ++j) {
            coords[j] = (i & (1 << j)) ? Tp(1) : Tp(0);
        }
        mesh.vertices.push_back(new Vertex<Tp, D>(coords, i));
    }

    std::vector<int> perm(D);
    for (int i = 0; i < D; ++i) perm[i] = i;

    Element<Tp, D> master = Element<Tp, D>::get_master_simplex();
    std::vector<Vertex<Tp, D>*> simplex_verts(D+1);
    do {
        simplex_verts[0] = mesh.vertices[0];  
        for (int i = 0; i < D; ++i) {
            int vertex_index = 0;
            for (int j = 0; j <= i; ++j) {
                vertex_index |= (1 << perm[j]);
            }
            assert(vertex_index >= 0 && vertex_index < mesh.vertices.size());
            simplex_verts[i + 1] = mesh.vertices[vertex_index];
        }
        Element<Tp, D> *simplex = new Element<Tp, D>(Element<Tp, D>::get_simplex(simplex_verts, master));
        mesh.elements.push_back(simplex);
    } while (std::next_permutation(perm.begin(), perm.end()));

    return mesh;
}
        

#endif /* MESH_HPP */