#include <iostream>
#include "matrix.h"

int main () {

    const uint32_t N = 929;
    const double a1 = 8.0;
    const double a2 = -1.0;
    const double a3 = -1.0;
    
    Matrix A(N, N);
    Matrix b(1, N);

    A.fillDiagonal(a1);
    A.fillDiagonal(a2, 1);
    A.fillDiagonal(a3, -1);
    A.fillDiagonal(a2, 2);
    A.fillDiagonal(a3, -2);

    for (int32_t i = 0; i < N; i++) {
        b(0, i) = sin(5*i);
    }

    Matrix result = jacobi(A, b);

    std::cout << result << std::endl;
}
