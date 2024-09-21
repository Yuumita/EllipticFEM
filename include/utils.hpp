#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cassert>


#include <eigen3/Eigen/Dense>
#include "types.h"

/// @brief Describes the linear function L(x) = coeffs^T * x + constant
/// @tparam Tp The scalar type
template <typename Tp>
class LinearFunction {
public:
    VectorX<Tp> coeffs;
    Tp constant;
    LinearFunction(const Tp constant_) : coeffs(0), constant(constant_) {}
    LinearFunction(const std::vector<Tp> coeffs_) : coeffs(coeffs_), constant(0) {}
    LinearFunction(const std::vector<Tp> coeffs_, const Tp constant_) : coeffs(coeffs_), constant(constant_) {}


    VectorX<Tp> gradient(VectorX<Tp> point) { return gradient(); }
    VectorX<Tp> gradient() { 
        return coeffs.head(coeffs.size() - 1); 
    }

    Tp operator()(const VectorX<Tp> &point) const {
        return evaluate(point);
    }

    Tp evaluate(const VectorX<Tp>& point) const {
        return constant + coeffs.dot(point);
    }

};

/// @brief Describes the affine transformation T(x) = Ax + b
/// @tparam Tp The scalar type
template <typename Tp>
class AffineTransformation {
public:
    MatrixX<Tp> A;
    VectorX<Tp> b;

    AffineTransformation();
    AffineTransformation(MatrixX<Tp> A_, VectorX<Tp> b_) : A(A_), b(b_) {
        assert(A.rows() == A.cols() && A.rows() == b.rows());
    }

    VectorX<Tp> operator()(const VectorX<Tp> &p) const {
        return A * p + b;
    }

    AffineTransformation inverse() const {
        MatrixX<Tp> A_inv = A.inverse();
        VectorX<Tp> bb = - A_inv * b;
        return AffineTransformation(A_inv, bb);
    }

};


/// @brief Compose an affine transformation with a linear function
/// @tparam Tp The scalar type
/// @return The function (lf; at)(x) = lf(at(x)).
template <typename Tp>
LinearFunction<Tp> compose_affine_linear(const AffineTransformation<Tp> &at, const LinearFunction<Tp> &lf) {
    return LinearFunction<Tp>(at.A.transpose() * lf.coeffs, lf.coeffs.transpose() * at.b + lf.constant);
}


template <typename Tp, int D>
class BilinearForm {
private: 
    std::function<MatrixX<Tp>(const VectorX<Tp>&)> Bf;
public:
    BilinearForm();
    BilinearForm(std::function<MatrixX<Tp>(const VectorX<Tp>&)> _Bf) : Bf(_Bf) {}
    Tp operator()(const LinearFunction<Tp> &f, const LinearFunction<Tp> &g, const Element<Tp, D> &element) const {
        Tp avg = Tp(0);
        for(std::pair<VectorX<Tp>, Tp> &p: element.get_quadrature_points()) {
            MatrixX<Tp> B = Bf ? Bf(p.first) : MatrixX<Tp>::Identity(D, D); 
            avg += p.second * (f.gradient(p.first).transpose() * B * g.gradient(p.first));
        }
        return avg * element.get_volume();
    }
    MatrixX<Tp> operator()(const VectorX<Tp> &x) const {
        return Bf(x);
    }
};

template <typename Tp, int D>
class LinearForm {
private: 
    std::function<Tp(const VectorX<Tp>&)> Lf;
public:
    LinearForm();
    LinearForm(std::function<Tp(const VectorX<Tp>&)> _Lf) : Lf(_Lf) {}
    Tp operator()(const LinearFunction<Tp> &f, const Element<Tp, D> &element) const {
        Tp avg = Tp(0);
        for(std::pair<VectorX<Tp>, Tp> &p: element.get_quadrature_points()) {
            avg += p.second * Lf(p.first) * f(p.first);
        }
        return avg * element.get_volume();
    }
    Tp operator()(const VectorX<Tp> &x) const {
        return Lf(x);
    }
};

template <typename Tp>
Tp factorial(size_t n) {
    static std::vector<Tp> _fact(1, Tp(1));
    while(_fact.size() <= n) _fact.push_back(_fact.back() * Tp(_fact.size()));
    return _fact[n];
}

#endif /* UTILS_HPP */