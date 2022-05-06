#include <iostream>
#include <cmath>
#include <tuple>

class Matrix {
private:
    uint32_t rows = 0;
    uint32_t cols = 0;
    double **mat = nullptr;

    void alloc();
    void dealloc();

public:
    Matrix() {}
    Matrix(uint32_t, uint32_t);
    Matrix(const Matrix&);
    Matrix(std::tuple<uint32_t, uint32_t>);
    Matrix& operator=(const Matrix&);

    ~Matrix();

    uint32_t getRows() const { return rows; }
    uint32_t getCols() const { return cols; }
    bool isVector() const { return rows == 1; }

    std::tuple<uint32_t, uint32_t> shape() const { return std::make_tuple(rows, cols); }

    void fillDiagonal(double value, int32_t k=0);
    double norm();
    void ones();

    inline double& operator()(uint32_t x, uint32_t y) const { return mat[x][y]; }

    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);
    Matrix& operator*=(double);
    Matrix& operator*=(const Matrix&);

    friend std::ostream& operator<< (std::ostream& os, const Matrix& matrix);
};

Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, const Matrix&);

const Matrix jacobi(const Matrix&, const Matrix&);
const Matrix gauss_seidel(const Matrix&, const Matrix&);
const Matrix lu(const Matrix& A, const Matrix& b);
