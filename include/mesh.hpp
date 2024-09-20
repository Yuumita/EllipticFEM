#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <eigen3/Eigen/Dense>
#include "types.h"
#include "vertex.hpp"
#include "element.hpp"


template <typename Tp, int D>
class Mesh {
public:
    static Mesh get_unit_cube_triangulation();
private:
    std::vector<Vertex<Tp, D>*>  vertices;
    std::vector<Element<Tp, D>*> elements;
};

template <typename Tp, int D>
Mesh<Tp, D> Mesh<Tp, D>::get_unit_cube_triangulation() {
    Mesh<Tp, D> mesh = new Mesh;

    VectorX<Tp> coords(D);
    for (int i = 0; i < (1 << D); ++i) {
        for (int j = 0; j < D; ++j) {
            coords[j] = (i & (1 << j)) ? Tp(1) : Tp(0);
        }
        mesh.vertices.push_back(new Vertex<Tp, D>(coords));
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
            simplex_verts[i + 1] = mesh.vertices[vertex_index];
        }
        Element<Tp, D> *simplex = new Element<Tp, D>(simplex_verts, master);
        mesh.elements.push_back(simplex);
    } while (std::next_permutation(perm.begin(), perm.end()));

    return mesh;
}
        

#endif /* MESH_HPP */