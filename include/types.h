#ifndef TYPES_H
#define TYPES_H

#include <eigen3/Eigen/Dense>
#include <iostream> // temporary inclusion

template<typename Tp>
using VectorX = Eigen::Matrix<Tp, Eigen::Dynamic, 1>;

template<typename Tp, int D>
using Vector = Eigen::Matrix<Tp, D, 1>;

template<typename Tp, int R, int C>
using Matrix = Eigen::Matrix<Tp, R, C>;

template<typename Tp>
using MatrixX = Eigen::Matrix<Tp, Eigen::Dynamic, Eigen::Dynamic>;

// #define DEBUG

#endif /* TYPES_H */