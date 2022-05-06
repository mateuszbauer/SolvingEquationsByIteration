#include <iostream>
#include <string.h>
#include <chrono>
#include "matrix.h"

#define JACOBI 0x1
#define GAUSS_SEIDEL 0x2
#define LU 0x4
#define ALL 0xff

int main (int argc, char* argv[]) {

    uint8_t opt = 0;

    for (size_t i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-jacobi") == 0) {
            opt |= JACOBI;
        }
        if (strcmp(argv[i], "-gs") == 0) {
            opt |= GAUSS_SEIDEL;
        }
        if (strcmp(argv[i], "-lu") == 0) {
            opt |= LU;
        }
    }

    if (opt == 0) {
        opt = ALL;
    }

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

    if (opt & GAUSS_SEIDEL) {
        std::cout << "Solving using Gauss-Seidel method...\n";
        auto start = std::chrono::steady_clock::now();
        Matrix result_gs = gauss_seidel(A, b);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << " ms\n\n";
    }

    if (opt & JACOBI) {
        std::cout << "Solving using Jacobi method...\n";
        auto start = std::chrono::steady_clock::now();
        Matrix result_jacobi = jacobi(A, b);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << " ms\n\n";
    }

    if (opt & LU) {
        std::cout << "Solving using LU factorization method...\n";
        auto start = std::chrono::steady_clock::now();
        Matrix result_lu = lu(A, b);
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << " ms\n\n";
    }

    return 0;
}
