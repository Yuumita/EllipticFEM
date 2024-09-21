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


template <typename Tp>
class BilinearForm {
public:
    MatrixX<Tp> B;
    BilinearForm();
    BilinearForm(MatrixX<Tp> B_) : B(B_) {}
    Tp operator()(const VectorX<Tp> &x, const VectorX<Tp> &y) const {
        return (x.transpose() * B) * y;
    }
};

template <typename Tp>
class LinearForm {
private: 
    VectorX<Tp> L_T;
public:
    LinearForm();
    LinearForm(const VectorX<Tp> &L_) : L_T(L_.transpose()) {}
    Tp operator()(const VectorX<Tp> &x) const {
        return L_T * x;
    }
};

template <typename Tp, int D>
Tp bilinear_quadrature_simplex(const BilinearForm<Tp> &B, const LinearFunction<Tp> &f, const LinearFunction<Tp> &g, const Element<Tp, D> &S) {
    VectorX<Tp> c = S.get_centroid();
    return S.get_volume() * B(f.gradient(c), g.gradient(c));
}

template <typename Tp, int D>
Tp linear_quadrature_simplex(const LinearForm<Tp> &L, const LinearFunction<Tp> &f, const Element<Tp, D> &S) {
    VectorX<Tp> c = S.get_centroid();
    return S.get_volume() * L(f(c));
}

template <typename Tp>
Tp factorial(size_t n) {
    static std::vector<Tp> _fact(1, Tp(1));
    while(_fact.size() <= n) _fact.push_back(_fact.back() * Tp(_fact.size()));
    return _fact[n];
}

#endif /* UTILS_HPP */