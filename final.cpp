#include <vector>
#include <iostream>

using namespace std;

// vector class comtable with doubles

template <typename T>
class Matrix
{
public:
    // initating
    Matrix();
    Matrix(int rows, int cols);
    Matrix(vector<vector<T>> &matrix);

    // operators
    T &operator()(int row, int col);
    const T &operator()(int row, int col) const;
    Matrix<T> &operator=(const Matrix<T> &other_matrix);
    friend ostream &operator<<(ostream &output, const Matrix<t> &matrix);

    void T();              // Transpose inplace;
    Matrix<T> transpose(); // Transpose returned;

private:
    int _rows;
    int _cols;

    vector<vector<T>> _matrix;
};

// initiating

template <typename T>
Matrix<T>::Matrix() : _rows(0), _cols(0)
{
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols) : _rows(rows), _cols(cols)
{
    vector<T> row_vector(cols);
    _matrix.resize(rows);

    // i = rows
    // j = cols

    for (int i = 0; i < rows; i++)
    {
        _matrix.push_back(row_vector);
    }
}

template <typename T>
Matrix<T>::Matrix(vector<vector<T>> matrix)
    _rows(matrix.size()),
    _cols(matrix[0].size()), _matrix(matrix) {}

/// operators

// set values
template <typename T>
T &Matrix<T>::operator()(int row, int col)
{
    if (row >= _rows || col >= _cols || row < 0 || col < 0)
    {
        throw out_of_range("Index out of bounds");
    }
    return _matrix[row][col];
}

// access (read) the component
template <typename T>
const T &Matrix<T>::operator()(int row, int col) const
{
    if (row >= _rows || col >= _cols || row < 0 || col < 0)
    {
        throw out_of_range("Index out of bounds");
    }
    return _matrix[row][col];
}

// set a defined matrix calss equal to another matrix.

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &other_matrix)
{
    _rows = other_matrix._rows;
    _cols = other_matrix._cols;
    _matrix = other_matrix._matrix;

    return *this;
}

// for printing
template <typename T>
ostream &operator<<(ostream &output, const Matrix<T> &matrix)
{
    for (int i = 0; i < matrix._rows; i++)
    {
        for (int j = 0; j < matrix._cols; ++j)
        {
            output << matrix._matrix[i][j] << " ";
        }
        output << endl;
    }
    return output;
}

/// Transpse functions

// returns transpose
template <typename T>
Matrix<T> Matrix<T>::transpose()
{
    Matrix<T> matrix_T(_cols, _rows);
    for (int j = 0; j < _cols; j++)
    {
        for (int i = 0; i < _rows; i++)
        {
            matrix_T._matrix[j][i] = _matrix[i][j];
        }
    }
    return matrix_T;
}

// sets curren tmatrix equal to its transpose
template <typename T>
void Matrix<T>::T()
{
    *this = transpose();
}

// Jacobian Method with pivoting

pair
