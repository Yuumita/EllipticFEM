#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cassert>


#include <eigen3/Eigen/Dense>
#include "types.h"

/// @brief Describes the linear function L(x) = coeffs^T * x + constant
/// @tparam Tp The scalar type
template <typename Tp, int D>
class LinearFunction {
public:
    Vector<Tp, D> coeffs;
    Tp constant;
    LinearFunction() : coeffs(Vector<Tp, D>::Zero()), constant(Tp(0)) {}
    LinearFunction(const Tp _constant) : coeffs(Vector<Tp, D>::Zero()), constant(_constant) {}
    LinearFunction(const Vector<Tp, D> _coeffs) : coeffs(_coeffs), constant(0) {}
    LinearFunction(const Vector<Tp, D> _coeffs, const Tp _constant) : coeffs(_coeffs), constant(_constant) {}


    Vector<Tp, D> gradient(Vector<Tp, D> point) const { return gradient(); }
    Vector<Tp, D> gradient() const {
        return coeffs;
    }

    Tp operator()(const Vector<Tp, D> &point) const {
        return constant + coeffs.dot(point);
    }

};

/// @brief Describes the affine transformation T(x) = Ax + b
/// @tparam Tp The scalar type
template <typename Tp, int D>
class AffineTransformation {
public:
    Matrix<Tp, D, D> A;
    Vector<Tp, D> b;

    AffineTransformation() : A(Matrix<Tp, D, D>::Identity()), b(Vector<Tp, D>::Zero()) {}
    AffineTransformation(Matrix<Tp, D, D> _A, Vector<Tp, D> _b) : A(_A), b(_b) { }

    Vector<Tp, D> operator()(const Vector<Tp, D> &p) const {
        return A * p + b;
    }

    AffineTransformation inverse() const {
        Matrix<Tp, D, D> A_inv = A.inverse();
        Vector<Tp, D> bb = - A_inv * b;
        AffineTransformation ret = AffineTransformation(A_inv, bb);
        return ret;
    }

};


/// @brief Compose an affine transformation with a linear function
/// @tparam Tp The scalar type
/// @return The function (lf; at)(x) = lf(at(x)).
template <typename Tp, int D>
LinearFunction<Tp, D> compose_affine_linear(const AffineTransformation<Tp, D> &at, const LinearFunction<Tp, D> &lf) {
    assert(at.A.cols() == lf.coeffs.size());
    Vector<Tp, D> new_coeffs = at.A.transpose() * lf.coeffs;
    Tp new_constant = lf.coeffs.dot(at.b) + lf.constant;
    return LinearFunction<Tp, D>(new_coeffs, new_constant);
}


template <typename Tp, int D>
class BilinearForm {
private: 
    std::function<Matrix<Tp, D, D>(const Vector<Tp, D>&)> Bf;
public:
    BilinearForm() : Bf([](const Vector<Tp, D>&) { return Matrix<Tp, D, D>::Identity(); }) {};
    BilinearForm(std::function<Matrix<Tp, D, D>(const Vector<Tp, D>&)> _Bf) : Bf(_Bf) {}
    Tp operator()(const LinearFunction<Tp, D> &f, const LinearFunction<Tp, D> &g, const Element<Tp, D> &element) const {
        Tp avg = Tp(0);
        for(std::pair<Vector<Tp, D>, Tp> &p: element.get_quadrature_points_1()) {
            Matrix<Tp, D, D> B = Bf ? Bf(p.first) : Matrix<Tp, D, D>::Identity(); 
            avg += p.second * Tp((f.gradient(p.first).transpose() * B) * g.gradient(p.first));
        }
        return avg;
    }
    Matrix<Tp, D, D> operator()(const Vector<Tp, D> &x) const { return Bf(x); }
};

template <typename Tp, int D>
class LinearForm {
private: 
    std::function<Tp(const Vector<Tp, D>&)> Lf;
public:
    LinearForm() : Lf([](const Vector<Tp, D>&) { return Tp(0); }) {}
    LinearForm(std::function<Tp(const Vector<Tp, D>&)> _Lf) : Lf(_Lf) {}
    Tp operator()(const LinearFunction<Tp, D> &f, const Element<Tp, D> &element) const {
        Tp avg = Tp(0); 
        for(std::pair<Vector<Tp, D>, Tp> &p: element.get_quadrature_points_rand()) {
            avg += p.second * Lf(p.first) * f(p.first);
        }
        return avg;
    }
    Tp operator()(const Vector<Tp, D> &x) const { return Lf(x); }
};

template <typename Tp>
Tp factorial(size_t n) {
    static std::vector<Tp> _fact(1, Tp(1));
    while(_fact.size() <= n) _fact.push_back(_fact.back() * Tp(_fact.size()));
    return _fact[n];
}

#endif /* UTILS_HPP */