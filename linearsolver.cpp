#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <random>
#include <ctime>
#include <chrono>
#include <algorithm>

std::vector<double> backwardsub(const std::vector<std::vector<double>> &U,
                                const std::vector<double> &b)
{
    int n = b.size();
    std::vector<double> x(n);

    x[n - 1] = b[n - 1] / U[n - 1][n - 1];

    for (int i = n - 2; i >= 0; --i)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j)
        {
            sum += U[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / U[i][i];
    }

    return x;
}

std::vector<double> forwardsub(const std::vector<std::vector<double>> &L,
                               const std::vector<double> &b)
{
    int n = b.size();
    std::vector<double> x(n);

    x[0] = b[0] / L[0][0];

    for (int i = 1; i < n; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
        {
            sum += L[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / L[i][i];
    }

    return x;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> gauselim(const std::vector<std::vector<double>> &A)
{
    int n = A.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U = A;

    for (int j = 0; j < n; ++j)
    {
        for (int i = j + 1; i < n; ++i)
        {
            if (U[j][j] == 0.0)
            {
                std::cout << "division by zero" << std::endl;
            }
            L[i][j] = U[i][j] / U[j][j];
            U[i][j] = 0.0;
            for (int k = j + 1; k < n; ++k)
            {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
    }

    for (int i = 0; i < n; ++i)
        L[i][i] = 1.0;

    return std::make_pair(U, L);
}

std::vector<std::vector<double>> Cholesky_Fact(const std::vector<std::vector<double>> &A)
{
    int n = A.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += L[i][k] * L[j][k];
            }
            if (i == j)
            {
                // Diagonal elements
                if (A[i][i] - sum < 0)
                {
                    throw std::runtime_error("Matrix is not positive definite.");
                }
                L[i][j] = sqrt(A[i][i] - sum);
            }
            else
            {
                // Off-diagonal elements
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }

    return L;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> GramSchmidt_Clasiscal(const std::vector<std::vector<double>> &A)
{
    int m = A.size();
    int n = A[0].size();

    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> Q(m, std::vector<double>(n, 0.0));

    std::vector<std::vector<double>> v(m, std::vector<double>(n, 0.0));

    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            v[i][j] = A[i][j];
        }

        for (int i = 0; i < j; ++i)
        {
            double dot_product = 0.0;
            for (int k = 0; k < m; ++k)
            {
                dot_product += Q[k][i] * v[k][j];
            }
            R[i][j] = dot_product;

            for (int k = 0; k < m; ++k)
            {
                v[k][j] -= R[i][j] * Q[k][i];
            }
        }

        double norm_vj = 0.0;
        for (int k = 0; k < m; ++k)
        {
            norm_vj += std::pow(v[k][j], 2);
        }
        R[j][j] = std::sqrt(norm_vj);

        for (int k = 0; k < m; ++k)
        {
            Q[k][j] = v[k][j] / R[j][j];
        }
    }

    return std::make_pair(Q, R);
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> GramSchmidt_Mod(const std::vector<std::vector<double>> &A)
{
    int m = A.size();
    int n = A[0].size();

    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> Q(m, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> V = A;

    for (int i = 0; i < n; ++i)
    {

        double norm = 0.0;
        for (int j = 0; j < m; ++j)
        {
            norm += V[j][i] * V[j][i];
        }
        R[i][i] = std::sqrt(norm);

        for (int j = 0; j < m; ++j)
        {
            Q[j][i] = V[j][i] / R[i][i];
        }

        for (int j = i + 1; j < n; ++j)
        {
            R[i][j] = 0.0;
            for (int k = 0; k < m; ++k)
            {
                R[i][j] += Q[k][i] * V[k][j];
            }
            for (int k = 0; k < m; ++k)
            {
                V[k][j] -= R[i][j] * Q[k][i];
            }
        }
    }

    return std::make_pair(Q, R);
}

double norm(const std::vector<double> &v)
{
    double sum = 0.0;
    for (double i : v)
    {
        sum += i * i;
    }
    return std::sqrt(sum);
}

void Householder(std::vector<std::vector<double>> &R, std::vector<double> &bb)
{
    int m = R.size();
    int n = R[0].size();

    for (int k = 0; k < n; ++k)
    {
        std::vector<double> x(m - k, 0);
        for (int i = k; i < m; ++i)
        {
            x[i - k] = R[i][k];
        }

        double x_norm = norm(x);
        double alpha = (x[0] >= 0) ? -x_norm : x_norm;
        x[0] -= alpha;

        double v_norm = norm(x);
        for (auto &vi : x)
        {
            vi /= v_norm; // normalizeing vector v ?
        }

        // appling transformation to all rows below and including k
        for (int i = k; i < m; ++i)
        {
            double dot = 0.0;
            for (int j = k; j < n; ++j)
            {
                dot += x[i - k] * R[i][j];
            }
            for (int j = k; j < n; ++j)
            {
                R[i][j] -= 2 * dot * x[i - k];
            }
        }

        // applying transformation to vector b
        double dot_b = 0.0;
        for (int i = k; i < m; ++i)
        {
            dot_b += x[i - k] * bb[i];
        }
        for (int i = k; i < m; ++i)
        {
            bb[i] -= 2 * dot_b * x[i - k];
        }
    }
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> partialpivot(const std::vector<std::vector<double>> &A)
{
    int n = A.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> P(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U = A;

    for (int i = 0; i < n; ++i)
    {
        P[i][i] = 1.0;
        L[i][i] = 1.0;
    }

    for (int j = 0; j < n - 1; ++j)
    {
        int max_index = j;
        double max_value = std::abs(U[j][j]);
        for (int i = j + 1; i < n; ++i)
        {
            if (std::abs(U[i][j]) > max_value)
            {
                max_index = i;
                max_value = std::abs(U[i][j]);
            }
        }
        if (max_index != j)
        {
            std::swap(U[j], U[max_index]);
            std::swap(P[j], P[max_index]);
            for (int k = 0; k < j; ++k)
            { // Correcting the lower part of L
                std::swap(L[j][k], L[max_index][k]);
            }
        }
        for (int i = j + 1; i < n; ++i)
        {
            L[i][j] = U[i][j] / U[j][j];
            U[i][j] = 0;
            for (int k = j + 1; k < n; ++k)
            {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
    }

    return std::make_tuple(U, L, P);
}

//// Functions for solving linear systems of equations Ax=b;

std::vector<double> linearsolver(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{

    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> UL = gauselim(A); // [U, L]

    std::vector<std::vector<double>> U = std::get<0>(UL);
    std::vector<std::vector<double>> L = std::get<1>(UL);

    std::vector<double> y = forwardsub(L, b);
    std::vector<double> x = backwardsub(U, y);

    return x;
}

std::vector<double> linearsolver_ppivot(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> ULP = partialpivot(A); // [U, L]

    std::vector<std::vector<double>> U = std::get<0>(ULP);
    std::vector<std::vector<double>> L = std::get<1>(ULP);
    std::vector<std::vector<double>> P = std::get<2>(ULP);

    std::vector<double> Pb(P.size());
    for (int i = 0; i < P.size(); i++)
    {
        double sum = 0;
        for (int j = 0; j < b.size(); j++)
        {
            sum += P[i][j] * b[j];
        }

        Pb[i] = sum;
    }

    std::vector<double> y = forwardsub(L, Pb);
    std::vector<double> x = backwardsub(U, y);

    return x;
}

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &matrix)
{
    int m = matrix.size();
    int n = matrix[0].size();
    std::vector<std::vector<double>> transposed(n, std::vector<double>(m, 0.0));
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            transposed[j][i] = matrix[i][j];
        }
    }
    return transposed;
}
std::vector<double> multiply(const std::vector<std::vector<double>> &matrix, const std::vector<double> &vec)
{
    int m = matrix.size();
    std::vector<double> result(m, 0.0);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < vec.size(); ++j)
        {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}

std::vector<double> GSC_LSQ_Solver(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{

    auto [Q, R] = GramSchmidt_Clasiscal(A);

    auto Q_transposed = transpose(Q);
    auto transformed_b = multiply(Q_transposed, b);

    auto X_lsq = backwardsub(R, transformed_b);

    return X_lsq;
}

std::vector<double> GSM_LSQ_Solver(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{

    auto [Q, R] = GramSchmidt_Mod(A);

    auto Q_transposed = transpose(Q);
    auto transformed_b = multiply(Q_transposed, b);

    auto X_lsq = backwardsub(R, transformed_b);

    return X_lsq;
}

std::vector<double> LSQsolver_Householder(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{
    std::vector<std::vector<double>> R = A;
    std::vector<double> bb = b;

    Householder(R, bb);
    return backwardsub(R, bb);
}

std::vector<double> LinearSolver_Chol(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{
    std::vector<std::vector<double>> RL = Cholesky_Fact(A);
    std::vector<std::vector<double>> RU = transpose(RL);
    std::vector<double> y = forwardsub(RL, b);
    std::vector<double> x = backwardsub(RU, y);
    return x;
}

/// Tests and Expirements
std::pair<std::vector<std::vector<double>>, std::vector<double>> generateRandomSystem(int size)
{
    std::vector<std::vector<double>> A(size, std::vector<double>(size));
    std::vector<double> b(size);
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            A[i][j] = dist(rng);
        }
        b[i] = dist(rng);
    }

    return {A, b};
}

// for positive sysmetrix matrics
std::pair<std::vector<std::vector<double>>, std::vector<double>> generate_random_PS(int size)
{
    std::vector<std::vector<double>> A(size, std::vector<double>(size));
    std::vector<double> b(size);
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (int i = 0; i < size; ++i)
    {
        for (int j = i; j < size; ++j)
        {
            A[i][j] = dist(rng);
            A[j][i] = A[i][j];
        }
        b[i] = dist(rng);
    }

    double max_value = 0;
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            max_value = std::max(max_value, std::abs(A[i][j]));
        }
    }
    for (int i = 0; i < size; ++i)
    {
        A[i][i] += max_value * size;
    }

    return {A, b};
}

int main(void)
{

    std::vector<int> sizes = {5, 10, 20, 50, 75, 100, 200, 500, 1000};
    // std::vector<int> sizes = {5, 10, 20, 50};

    for (int size : sizes)
    {
        auto [A, b] = generate_random_PS(size);
        using namespace std::chrono;
        steady_clock::time_point start;
        steady_clock::time_point end;

        std::vector<double> x;

        // linearsolver
        start = steady_clock::now();
        x = linearsolver(A, b);
        end = steady_clock::now();
        std::cout << "Time for linearsolver with size " << size << ": " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

        // linearsolver_ppivot
        start = steady_clock::now();
        x = linearsolver_ppivot(A, b);
        end = steady_clock::now();
        std::cout << "Time for linearsolver_ppivot with size " << size << ": " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

        // GSC_LSQ_Solver
        start = steady_clock::now();
        x = GSC_LSQ_Solver(A, b);
        end = steady_clock::now();
        std::cout << "Time for GSC_LSQ_Solver with size " << size << ": " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

        // GSM_LSQ_Solver
        start = steady_clock::now();
        x = GSM_LSQ_Solver(A, b);
        end = steady_clock::now();
        std::cout << "Time for GSM_LSQ_Solver with size " << size << ": " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

        // LSQsolver_Householder
        start = steady_clock::now();
        x = LSQsolver_Householder(A, b);
        end = steady_clock::now();
        std::cout << "Time for LSQsolver_Householder with size " << size << ": " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

        start = steady_clock::now();
        x = LinearSolver_Chol(A, b);
        end = steady_clock::now();
        std::cout << "Time for LinearSolver_Chol with size " << size << ": " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

        std::cout << "\n";
    }

    auto [A, b] = generate_random_PS(10);

    std::vector<double> x = LinearSolver_Chol(A, b);

    std::cout << "Solution for linearsolver at size 10:\n";
    x = linearsolver(A, b);
    for (double val : x)
    {
        std::cout << val << " ";
    }
    std::cout << "\n";

    x = linearsolver_ppivot(A, b);
    std::cout << "Solution for linearsolver_ppivot at size 10:\n";
    for (double val : x)
    {
        std::cout << val << " ";
    }
    std::cout << "\n";

    x = GSC_LSQ_Solver(A, b);
    std::cout << "Solution for GSC_LSQ_Solver at size 10:\n";
    for (double val : x)
    {
        std::cout << val << " ";
    }
    std::cout << "\n";

    x = GSM_LSQ_Solver(A, b);
    std::cout << "Solution for GSM_LSQ_Solver at size 10:\n";
    for (double val : x)
    {
        std::cout << val << " ";
    }
    std::cout << "\n";

    // x = LSQsolver_Householder(A, b);
    // std::cout << "Solution for LSQsolver_Householder at size 10:\n";
    // for (double val : x)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << "\n";

    x = LinearSolver_Chol(A, b);
    std::cout << "Solution for LinearSolver_Chol at size 10:\n";
    for (double val : x)
    {
        std::cout << val << " ";
    }
    std::cout << "\n";

    return 0;
}
