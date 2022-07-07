#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>

#define _USE_MATH_DEFINES

int main()
{
    // Basic Example of cpp
    std::cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    std::cout << a << std::endl;
    std::cout << a / b << std::endl;
    std::cout << std::sqrt(b) << std::endl;
    std::cout << std::acos(-1) << std::endl;
    std::cout << std::sin(30.0 / 180.0 * acos(-1)) << std::endl;
    // Example of vector
    std::cout << "Example of vector \n";
    // vector definition
    Eigen::Vector3f v(1.0f, 2.0f, 3.0f);
    Eigen::Vector3f w(1.0f, 0.0f, 0.0f);
    // vector output
    std::cout << "Example of output \n";
    std::cout << v << std::endl;
    // vector add
    std::cout << "Example of add \n";
    std::cout << v + w << std::endl;
    // vector scalar multiply
    std::cout << "Example of scalar multiply \n";
    std::cout << v * 3.0f << std::endl;
    std::cout << 2.0f * v << std::endl;
    // vector dot product
    std::cout << "Example of vector dot product \n";
    std::cout << v.transpose() * w << std::endl;
    // Example of matrix
    std::cout << "Example of matrix \n";
    // matrix definition
    Eigen::Matrix3f i, j;
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    std::cout << "Example of output \n";
    std::cout << i << std::endl;
    // matrix add i + j
    std::cout << "Example of matrix add \n";
    std::cout << i + j << std::endl;
    // matrix scalar multiply i * 2.0
    std::cout << "Example of matrix scalar multiply \n";
    std::cout << i * 2.0f << std::endl;
    // matrix multiply i * j
    std::cout << "Example of matrix multiply \n";
    std::cout << i * j << std::endl;
    // matrix multiply vector i * v
    std::cout << "Example of matrix multiply vector \n";
    std::cout << i * v << std::endl;
    // homework
    Eigen::Vector3d p(2.0f, 1.0f, 1.0f);
    Eigen::Matrix3d rotate;
    Eigen::Matrix3d trans;
    double theta = 45.0 / 180.0 * M_PI;
    rotate << cos(theta), -1.0 * sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
    trans << 1, 0, 1, 0, 1, 2, 0, 0, 1;
    p = trans * rotate * p;
    std::cout << p << std::endl;
    return 0;
}