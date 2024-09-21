#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <vector>
#include "types.h"


template <typename Tp, int D>
class Element;

/// @brief 
/// @tparam Tp 
/// @tparam D 
template <typename Tp, int D>
class Vertex {
private:
    int index;
    VectorX<Tp> coords;
    // std::vector<Element<Tp, D>*> elements;
public:
    Vertex();
    Vertex(const Tp &v0, const int _index): index(_index) {
        coords.resize(D);
        for(int i = 0; i < D; i++) coords[i] = v0;
    }
    Vertex(const VectorX<Tp> &p): coords(p) {}
    Vertex(const VectorX<Tp> &p, const int _index): coords(p), index(_index) {}

    int get_index()  { return index; }
    VectorX<Tp> get_coords() { return coords; }
};



#endif /* VERTEX_HPP */