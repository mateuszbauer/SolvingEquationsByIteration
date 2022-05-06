#include <iostream>
#include <cassert>
#include "matrix.h"

Matrix::Matrix(uint32_t rows, uint32_t cols) : rows(rows), cols(cols) {
    alloc();
}

Matrix::Matrix(const Matrix& m) : rows(m.rows), cols(m.cols) {
    alloc();

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat[i][j] = m.mat[i][j];
        }
    }
}

Matrix::Matrix(std::tuple<uint32_t, uint32_t> shape) : Matrix(std::get<0>(shape), std::get<1>(shape)) {}

void Matrix::alloc () {
    assert (mat == nullptr);
    mat = new double*[rows];

    for (size_t i = 0; i < rows; i++) {
        mat[i] = new double[cols];
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat[i][j] = 0;
        }
    }
}

Matrix& Matrix::operator=(const Matrix& m) {
    if (this == &m) {
        return *this;
    }
    if (this->shape() != m.shape()) {
        if (mat != nullptr) {
            dealloc();
        }
        rows = m.rows;
        cols = m.cols;
        alloc();
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat[i][j] = m.mat[i][j];
        }
    }

    return *this;
}
    

void Matrix::dealloc () {
    assert (mat != nullptr);
    for (size_t i = 0; i < rows; i++) {
        delete[] mat[i];
    }

    delete[] mat;
    mat = nullptr;
}

Matrix::~Matrix() {
    dealloc();
}

Matrix& Matrix::operator+=(const Matrix& m) {
    assert (this->shape() == m.shape());
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat[i][j] += m.mat[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& m) {
    assert (this->shape() == m.shape());
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat[i][j] -= m.mat[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(double n) {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat[i][j] *= n;
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& m) {
    if (this->getCols() != m.getRows()) {
        assert(m.isVector());
        assert(this->getCols() == m.getCols());
    }

    Matrix tmp(*this);
    dealloc();

    if (m.isVector()) {
        rows = 1;
        alloc();

        for (size_t i = 0; i < tmp.cols; i++) {
            for (size_t j = 0; j < tmp.rows; j++) {
                mat[0][i] += tmp(j, i) * m(0, j);
            }
        }
        return *this;
    }
    else {
        cols = m.cols;
        alloc();

        for (size_t i = 0; i < tmp.rows; i++) {
            for (size_t j = 0; j < m.cols; j++) {

                for (size_t k = 0; k < m.rows; k++) {
                    mat[i][j] += tmp(i, k) * m(k, j);
                }
            }
        }
        return *this;
    }
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    os << "Matrix shape = (" << matrix.rows << ", " << matrix.cols << ")\n";
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.cols; j++) {
            os << matrix.mat[i][j] << " ";
        }
        os << std::endl;
    }

    return os;
}

Matrix operator+(const Matrix& m1, const Matrix& m2) {
    assert (m1.shape() == m2.shape());
    Matrix tmp(m1);
    return tmp += m2;
}

Matrix operator-(const Matrix& m1, const Matrix& m2) {
    assert(m1.shape() == m2.shape());
    Matrix tmp(m1);
    return tmp -= m2;
}

Matrix operator*(const Matrix& m1, const Matrix& m2) {
    if (m1.getCols() != m2.getRows()) {
        assert(m2.isVector());
        assert(m1.getCols() == m2.getCols());
    }

    Matrix tmp(m1);
    return tmp *= m2;
}

void Matrix::fillDiagonal(double value, int32_t k) {
    int32_t i = 0;
    int32_t j = 0;

    (k >= 0) ? j += k : i -=k;

    while (i < rows && j < cols) {
        mat[i][j] = value;
        i++;
        j++;
    }
}

double Matrix::norm() {
    double ret = 0;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            ret += pow(mat[i][j], 2);
        }
    }

    return sqrt(ret);
}

void Matrix::ones() {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat[i][j] = 1.0;
        }
    }
}

const Matrix jacobi(const Matrix& A, const Matrix& b) {
    assert(A.getCols() == A.getRows());
    assert(b.isVector());
    assert(A.getCols() == b.getCols());

    Matrix x(b.shape());
    Matrix x_new(b.shape());
    Matrix residuum(b.shape());

    do {
        for (size_t i = 0; i < A.getRows(); i++) {
            double k = 0;
            for (size_t j = 0; j < A.getCols(); j++) {
                if (i != j) {
                    k += A(i, j) * x(0, j);
                }
            }
            x_new(0, i) = (b(0, i) - k) / A(i, i);
        }
        x = x_new;
        residuum = (A * x) - b;       
    } while (residuum.norm() > 1e-9);

    return x;
}

const Matrix gauss_seidel(const Matrix& A, const Matrix& b) {
    assert(A.getCols() == A.getRows());
    assert(b.isVector());
    assert(A.getCols() == b.getCols());

    Matrix x(b.shape());
    Matrix x_new(b.shape());
    Matrix residuum(b.shape());

    do {
        for (size_t i = 0; i < A.getRows(); i++) {
            double k = 0;
            for (size_t j = 0; j < i; j++) {
                k += A(i, j) * x_new(0, j);
            }
            for (size_t j = i + 1; j < A.getCols(); j++) {
                k += A(i, j) * x(0, j);
            }
            x_new(0, i) = (b(0, i) - k) / A(i, i);
        }
        x = x_new;
        residuum = (A * x) - b;
    } while (residuum.norm() > 1e-9);

    return x;
}
