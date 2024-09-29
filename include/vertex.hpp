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
    Vector<Tp, D> coords;
    // std::vector<Element<Tp, D>*> elements;
public:
    Vertex();
    Vertex(const Tp &v0) {
        coords.resize(D);
        for(int i = 0; i < D; i++) coords[i] = v0;
    }
    Vertex(const Tp &v0, const int _index): index(_index) {
        coords.resize(D);
        for(int i = 0; i < D; i++) coords[i] = v0;
    }
    Vertex(const Vector<Tp, D> &p): coords(p) {}
    Vertex(const Vector<Tp, D> &p, const int _index): coords(p), index(_index) {}

    int get_index() const { return index; }
    Vector<Tp, D>& get_coords() { return coords; }

    Tp &coefRef(size_t i)    { return coords[i]; }
    Tp &operator[](size_t i) { return coords[i]; }

    friend std::ostream& operator<<(std::ostream& os, const Vertex<Tp, D>& vertex) {
        os << "([";
        for (int i = 0; i < D; ++i) {
            os << vertex.coords[i];
            if (i < D - 1) os << ", ";
        }
        os << "], " << vertex.index << ")";
        return os;
    }

};



#endif /* VERTEX_HPP */