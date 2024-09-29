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
    static Mesh<Tp, D> get_unit_cube_triangulation();
    static Mesh<Tp, D> get_2d_unit_cube_triangulation(Tp h = 1);
    static Mesh<Tp, D> get_3d_unit_cube_triangulation(Tp h = 1);

    Vertex<Tp, D> *get_vertex(size_t i)   { return vertices[i]; }
    std::vector<Vertex<Tp, D>*> &get_vertices()   { return vertices; }
    size_t get_vertices_size()   { return vertices.size(); }

    Element<Tp, D> *get_element(size_t i) { return elements[i]; }
    std::vector<Element<Tp, D>*> &get_elements()   { return elements; }
    size_t get_elements_size()   { return elements.size(); }

};


template <typename Tp, int D>
Mesh<Tp, D> Mesh<Tp, D>::get_2d_unit_cube_triangulation(Tp h) {
    assert(D == 2);
    Mesh<Tp, 2> mesh;

    int k = static_cast<int>(1 / h);

    for(int x = 0; x <= k; x++) {
        for(int y = 0; y <= k; y++) {
            Vector<Tp, 2> coords;
            coords[0] = x * h;
            coords[1] = y * h;
            mesh.vertices.push_back(new Vertex<Tp, 2>(coords, x * (k+1) + y));
#ifdef DEBUG
            std::cerr << "Created vertex " << *mesh.vertices.back() << std::endl;
#endif
        }
    }


    Element<Tp, 2> master = Element<Tp, 2>::get_master_simplex();
    
    for(int x = 0; x < k; x++) {
        for(int y = 0; y < k; y++) {
            int base = x * (k+1) + y;
            std::vector<int> subsquare_indices = { base, base + 1, base + (k+1), base + (k+1) + 1 };

            std::vector<Vertex<Tp, 2>*> subsqr(4);
            for(int i = 0; i < 4; i++) {
                subsqr[i] = mesh.vertices[subsquare_indices[i]];
            }

            Vector<Tp, 2> coords;
            coords[0] = (subsqr[0]->get_coords()[0] + subsqr[3]->get_coords()[0]) / Tp(2);
            coords[1] = (subsqr[0]->get_coords()[1] + subsqr[3]->get_coords()[1]) / Tp(2);
            Vertex<Tp, 2> *c = new Vertex<Tp, 2>(coords, mesh.vertices.size());
            mesh.vertices.push_back(c);
#ifdef DEBUG
            std::cerr << "Created vertex: " << *c << std::endl;
#endif


            std::vector<std::vector<Vertex<Tp, 2>*>> simplices_verts = 
                {
                    {c, subsqr[0], subsqr[1]},
                    {c, subsqr[0], subsqr[2]},
                    {c, subsqr[1], subsqr[3]},
                    {c, subsqr[2], subsqr[3]}
                };

            for(int i = 0; i < 4; i++) {
                Element<Tp, 2> *simplex = new Element<Tp, 2>(Element<Tp, 2>::get_simplex(simplices_verts[i], master));
                mesh.elements.push_back(simplex);
#ifdef DEBUG
            std::cerr << "Created triangle between " << 
                "(" << simplices_verts[i][0]->get_coords()[0] << ", " << simplices_verts[i][0]->get_coords()[1] << ") - " << 
                "(" << simplices_verts[i][1]->get_coords()[0] << ", " << simplices_verts[i][1]->get_coords()[1] << ") - " << 
                "(" << simplices_verts[i][2]->get_coords()[0] << ", " << simplices_verts[i][2]->get_coords()[1] << ")" <<  std::endl;
#endif
            }
        }
    }

    return mesh;
}

template <typename Tp, int D>
Mesh<Tp, D> Mesh<Tp, D>::get_3d_unit_cube_triangulation(Tp h) {
    assert(D == 3);
    Mesh<Tp, 3> mesh;

    int k = static_cast<int>(1 / h);

    for(int x = 0; x <= k; x++) {
        for(int y = 0; y <= k; y++) {
            for(int z = 0; z <= k; z++) {
                Vector<Tp, 3> coords;
                coords[0] = x * h;
                coords[1] = y * h;
                coords[2] = z * h;
                mesh.vertices.push_back(new Vertex<Tp, 3>(coords, x * (k+1) * (k+1) + y * (k+1) + z));

            }
        }
    }


    Element<Tp, 3> master = Element<Tp, 3>::get_master_simplex();

    for(int x = 0; x < k; x++) {
        for(int y = 0; y < k; y++) {
            for(int z = 0; z < k; z++) {
                int base = x * (k+1) * (k+1) + y * (k+1) + z;
                std::vector<int> subcube_indices = {
                    base, base + 1, base + (k+1), base + (k+1) + 1, 
                    base + (k+1)*(k+1), base + (k+1)*(k+1) + 1, 
                    base + (k+1)*(k+1) + (k+1), base + (k+1)*(k+1) + (k+1) + 1
                };

                std::vector<Vertex<Tp, 3>*> subcube(8);
                for(int i = 0; i < 8; i++) {
                    subcube[i] = mesh.vertices[subcube_indices[i]];
                }

                std::vector<std::vector<Vertex<Tp, 3>*>> simplices_verts = 
                {
                    {subcube[0], subcube[1], subcube[1+2], subcube[1+2+4]},
                    {subcube[0], subcube[2], subcube[1+2], subcube[1+2+4]},

                    {subcube[0], subcube[1], subcube[1+4], subcube[1+2+4]},
                    {subcube[0], subcube[2], subcube[2+4], subcube[1+2+4]},

                    {subcube[0], subcube[4], subcube[1+4], subcube[1+2+4]},
                    {subcube[0], subcube[4], subcube[2+4], subcube[1+2+4]}
                };
                for(int i = 0; i < 6; i++) {
                    // std::vector<Vertex<Tp, 3>*> simplex_verts(3+1);
                    // for(int j = 0; j < 3+1; j++) {
                    //     simplex_verts[j] = simplices_verts[i][j];
                    // }
                    Element<Tp, 3> *simplex = new Element<Tp, 3>(Element<Tp, 3>::get_simplex(simplices_verts[i], master));
                    mesh.elements.push_back(simplex);
                }
            }
        }
    }

    return mesh;
}

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