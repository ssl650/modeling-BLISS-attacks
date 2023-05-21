#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main() {
    Matrix<long,  1, 10> vec;
    Matrix<long, 10, 10> mat;

    vec << 1, 5, 17, -3, -2, 6, 12, 4, 9, 26;

    mat = vec.transpose() * vec;

    Matrix<double, 10, 10> md = mat.template cast<double>();

    std::cout << "The matrix M = v^T * v is:\n" << md << std::endl;

    SelfAdjointEigenSolver<Matrix<double, 10, 10> > es;
    es.compute(md);

    std::cout << "The eigenvalues of M are:\n" << es.eigenvalues()
	      << std::endl;
    std::cout << "Eigenvector corresponding to the nonzero ev of M:\n"
	      << es.eigenvectors().col(9) << std::endl;

    std::cout << "Check:\n";
    Matrix<double, 10, 1> ev = es.eigenvectors().col(9), mev;
    mev = md * ev;
    std::cout << "norm(ev)   = " << ev.norm() << std::endl;
    std::cout << "norm(M*ev) = " << mev.norm() << std::endl;


    return 0;
}



