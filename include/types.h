#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Dense>

template<typename Tp>
using VectorX = Eigen::Matrix<Tp, Eigen::Dynamic, 1>;

template<typename Tp>
using MatrixX = Eigen::Matrix<Tp, Eigen::Dynamic, Eigen::Dynamic>;

#endif /* TYPES_H */