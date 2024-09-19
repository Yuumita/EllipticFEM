#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cassert>


#include <Eigen/Dense>
#include "types.h"

/// @brief Describes the linear function L(x) = coeffs^T * x + constant
/// @tparam Tp The scalar type
template <typename Tp>
class LinearFunction {
private:
    VectorX<Tp> coeffs;
    Tp constant;
public:
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
private:
    MatrixX<Tp> A;
    VectorX<Tp> b;
public:
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
    return LinearFunction<Tp>(at.A.T * lf.coeffs, lf.coeffs.T * at.b + lf.constant);
}


#endif /* UTILS_HPP */