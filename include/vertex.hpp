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
    // size_t index;
    // std::vector<Element<Tp, D>*> elements;
    VectorX<Tp> coords;
public:
    Vertex();
    Vertex(const Tp &v0) {
        coords.resize(D);
        for(int i = 0; i < D; i++) coords[i] = v0;
    }
    Vertex(const VectorX<Tp> &p);

    VectorX<Tp> get_coords() { return coords; }
};



#endif /* VERTEX_HPP */