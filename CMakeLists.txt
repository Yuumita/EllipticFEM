cmake_minimum_required(VERSION 3.10)
project(EllipticFEM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(EllipticFEM src/main.cpp)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(EllipticFEM Eigen3::Eigen)
