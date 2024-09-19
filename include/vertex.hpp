#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <vector>
#include "element.hpp"

#include "types.h"


template <typename Tp, int D>
class Element;

/// @brief 
/// @tparam Tp 
/// @tparam TPoint 
/// @tparam D 
template <typename Tp, int D>
class Vertex {
private:
    size_t index;
    std::vector<Element<Tp, D>*> elements;
    VectorX<Tp> coords;
public:
    Vertex();
    Vertex(const Tp &v0) {
        coords.resize(D);
        for(int i = 0; i < D; i++) coords[i] = v0;
    }
    Vertex(const VectorX<Tp> &p);
};



#endif /* VERTEX_HPP */