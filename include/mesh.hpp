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
    std::vector<Vertex<Tp, D>*>  inner_vertices;
    std::vector<Vertex<Tp, D>*>  boundary_vertices;
    std::vector<Element<Tp, D>*> elements;

public:
    static Mesh<Tp, D> get_unit_cube_triangulation(int divs = 1);
    static Mesh<Tp, D> get_2d_unit_cube_triangulation(int divs = 1);
    static Mesh<Tp, D> get_3d_unit_cube_triangulation(int divs = 1);


    std::vector<Vertex<Tp, D>*> &get_boundary_vertices()   { return boundary_vertices; }
    std::vector<Vertex<Tp, D>*> &get_inner_vertices()   { return inner_vertices; }
    size_t get_inner_vertices_size() const { return inner_vertices.size(); }
    size_t get_boundary_vertices_size() const { return inner_vertices.size(); }

    Vertex<Tp, D> *get_inner_vertex(size_t i)   { 
        if (i >= inner_vertices.size()) throw std::out_of_range("Inner vertex index out of range");
        return inner_vertices[i]; 
    }

    Vertex<Tp, D> *get_boundary_vertex(size_t i)   { 
        if (i >= boundary_vertices.size()) throw std::out_of_range("Boundary vertex index out of range");
        return boundary_vertices[i]; 
    }


    Element<Tp, D> *get_element(size_t i) { 
        if (i >= elements.size()) throw std::out_of_range("Element index out of range");
        return elements[i];
    }
    std::vector<Element<Tp, D>*> &get_elements()   { return elements; }
    size_t get_elements_size()   { return elements.size(); }

};


template <typename Tp, int D>
Mesh<Tp, D> Mesh<Tp, D>::get_2d_unit_cube_triangulation(int divs) {
    static_assert(D == 2, "get_2d_unit_cube_triangulation can only be called for 2D meshes.");

    Mesh<Tp, 2> mesh;

    Tp h = Tp(1) / static_cast<Tp>(divs);
    std::vector<Vertex<Tp, D>*> vertices;

    int index_c = 0;
    for(int x = 0; x <= divs; x++) {
        for(int y = 0; y <= divs; y++) {
            Vector<Tp, 2> coords;
            coords[0] = x * h;
            coords[1] = y * h;
            bool at_boundary = coords[0] == Tp(0) || coords[0] == Tp(1) || coords[1] == Tp(0) || coords[1] == Tp(1);
            Vertex<Tp, 2> *c = new Vertex<Tp, 2>(coords, (at_boundary ? -1 : index_c++));
            vertices.push_back(c);
            if(!at_boundary) mesh.inner_vertices.push_back(c);
            else mesh.boundary_vertices.push_back(c);
        }
    }


    Element<Tp, 2> master = Element<Tp, 2>::get_master_simplex();
    
    for(int x = 0; x < divs; x++) {
        for(int y = 0; y < divs; y++) {
            int base = x * (divs+1) + y;
            std::vector<int> subsquare_indices = { base, base + 1, base + (divs+1), base + (divs+1) + 1 };

            std::vector<Vertex<Tp, 2>*> subsqr(4);
            for(int i = 0; i < 4; i++) {
                subsqr[i] = vertices[subsquare_indices[i]];
            }

            Vector<Tp, 2> coords;
            coords[0] = (subsqr[0]->get_coords()[0] + subsqr[3]->get_coords()[0]) / Tp(2);
            coords[1] = (subsqr[0]->get_coords()[1] + subsqr[3]->get_coords()[1]) / Tp(2);
            bool at_boundary = coords[0] == Tp(0) || coords[0] == Tp(1) || coords[1] == Tp(0) || coords[1] == Tp(1);
            Vertex<Tp, 2> *c = new Vertex<Tp, 2>(coords, (!at_boundary ? index_c++ : -1));
            vertices.push_back(c);
            if(!at_boundary) mesh.inner_vertices.push_back(c);
            else mesh.boundary_vertices.push_back(c);


            std::vector<std::vector<Vertex<Tp, 2>*>> simplices_verts = 
                {
                    {c, subsqr[0], subsqr[1]},
                    {c, subsqr[0], subsqr[2]},
                    {c, subsqr[1], subsqr[3]},
                    {c, subsqr[2], subsqr[3]}
                };

            for(int i = 0; i < 4; i++) {
                Element<Tp, 2> *simplex = new Element<Tp, 2>(Element<Tp, 2>::get_simplex(simplices_verts[i], master));
                bool at_boundary = false;
                mesh.elements.push_back(simplex);
            }
        }
    }

    return std::move(mesh);
}

template <typename Tp, int D>
Mesh<Tp, D> Mesh<Tp, D>::get_unit_cube_triangulation(int divs) {
    Mesh<Tp, D> mesh;

    Tp h = Tp(1) / static_cast<Tp>(divs);
    std::vector<int> coord_indices(D, 0);

    int verts_size = std::pow(divs + 1, D);
    std::vector<Vertex<Tp, D>*> vertices(verts_size);
    int inner_verts_index = 0;

    for(int i = 0; i < verts_size; i++) {
        Vector<Tp, D> coords(D);
        int idx = i;
        bool at_boundary = false;
        for(int j = 0; j < D; j++) {
            coord_indices[j] = idx % (divs + 1);
            coords[j] = coord_indices[j] * h;
            if(coords[j] == Tp(0) || coords[j] == Tp(1))
                at_boundary = true;
            idx /= (divs + 1);
        }
        Vertex<Tp, D> *vert = new Vertex<Tp, D>(coords, (at_boundary ? -1 : inner_verts_index++));
        vertices[i] = vert;
        if(!at_boundary) mesh.inner_vertices.push_back(vert);
        else mesh.boundary_vertices.push_back(vert);
    }


    Element<Tp, D> master = Element<Tp, D>::get_master_simplex();

    #pragma omp parallel for
    for(int base_idx = 0; base_idx < verts_size; base_idx++) {
        bool valid_base = true;
        for(int j = 0; valid_base && j < D; j++)
            if(vertices[base_idx]->coefRef(j) == Tp(1)) valid_base = false;
        if(!valid_base) continue;

        Vector<int, D> base_indices;
        int idx = base_idx;
        for(int j = 0; j < D; j++) {
            base_indices[j] = idx % (divs + 1);
            idx /= (divs + 1);
        }

        std::vector<Vertex<Tp, D>*> cube_verts(1 << D);
        for(int i = 0; i < (1 << D); i++) {
            Vector<int, D> vertex_indices = base_indices;
            for(int j = 0; j < D; j++) {
                if(i & (1 << j)) vertex_indices[j]++;
            }

            int idx = 0;
            for(int j = D-1; j >= 0; j--) {
                idx = (divs + 1) * idx + vertex_indices[j];
            }
            cube_verts[i] = vertices[idx];
        }

        std::vector<int> perm(D);
        for (int i = 0; i < D; ++i) perm[i] = i;

        std::vector<Vertex<Tp, D>*> simplex_verts(D+1);
        do {
            simplex_verts[0] = cube_verts[0];  
            for (int i = 0; i < D; ++i) {
                int vertex_index = 0;
                for (int j = 0; j <= i; ++j) {
                    vertex_index |= (1 << perm[j]);
                }
                simplex_verts[i + 1] = cube_verts[vertex_index];
            }
            Element<Tp, D> *simplex = new Element<Tp, D>(Element<Tp, D>::get_simplex(simplex_verts, master));
            mesh.elements.push_back(simplex);
        } while (std::next_permutation(perm.begin(), perm.end()));

    }

    return std::move(mesh);
}

// template <typename Tp, int D>
// Mesh<Tp, D> Mesh<Tp, D>::get_unit_cube_triangulation() {
//     Mesh<Tp, D> mesh;
// 
//     Vector<Tp, D> coords(D);
//     for (int i = 0; i < (1 << D); ++i) {
//         for (int j = 0; j < D; ++j) {
//             coords[j] = (i & (1 << j)) ? Tp(1) : Tp(0);
//         }
//         mesh.vertices.push_back(new Vertex<Tp, D>(coords, i));
//     }
// 
//     std::vector<int> perm(D);
//     for (int i = 0; i < D; ++i) perm[i] = i;
// 
//     Element<Tp, D> master = Element<Tp, D>::get_master_simplex();
//     std::vector<Vertex<Tp, D>*> simplex_verts(D+1);
//     do {
//         simplex_verts[0] = mesh.vertices[0];  
//         for (int i = 0; i < D; ++i) {
//             int vertex_index = 0;
//             for (int j = 0; j <= i; ++j) {
//                 vertex_index |= (1 << perm[j]);
//             }
//             assert(vertex_index >= 0 && vertex_index < mesh.vertices.size());
//             simplex_verts[i + 1] = mesh.vertices[vertex_index];
//         }
//         Element<Tp, D> *simplex = new Element<Tp, D>(Element<Tp, D>::get_simplex(simplex_verts, master));
//         mesh.elements.push_back(simplex);
//     } while (std::next_permutation(perm.begin(), perm.end()));
// 
//     return std::move(mesh);
// }
        

#endif /* MESH_HPP */