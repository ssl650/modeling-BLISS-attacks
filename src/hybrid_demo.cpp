#include "Parameters.h"
#include "Setup.h"
#include "Sampler.h"
#include "KeyGen.h"
#include <iostream>
#include <NTL/LLL.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <map>
#include <random>
#include <algorithm>

inline double static calculateHermitFactor(double beta)
{
    return pow(beta * pow(M_PI * beta, 1 / beta) / (2 * M_PI * M_E),
               1 / (2 * (beta - 1)));
}

//TODO: reorganize parameters set
#define beta 10 // block size for BKZ
#define r N // meet-in-the-middle dimension
#define delta 0.75

int attackTime;

// Function to permute an NTL vector
void permuteVector(NTL::vec_ZZ& vector)
{
    // Perform permutation
    std::random_device rd;
    std::mt19937 gen(rd());
    shuffle(vector.begin(), vector.end(), gen);
}

void vectors_for_A(long y, NTL::vec_ZZ x, NTL::mat_ZZ& res) 
{
    // инициализация граничных значений
    double xp = static_cast<double>(y) / 2 - 1;
    double xn = static_cast<double>(y) / 2;
    long border_c = static_cast<long>(ceil(xp));
    long border_f = static_cast<long>(-floor(xn));

    // подсчет столбцов и определение постоянных значений
    NTL::vec_ZZ result;
    result.SetLength(N);
    long rows = 1;
    for (int i = 0; i <N ; i++)
    {
        if (x[i] <= border_c and x[i] >= border_f)
        {
            rows *= 2;
            result[i] = -1;
        }
        else if (x[i] > border_c)
            result[i] = 1;
        else if (x[i] < border_f)
            result[i] = 0;
    }
    res.SetDims(rows, N);

    // инициализация всех столбцов
    for (int i = 0; i < rows; i++)
    {
        int ind = 0;
        for (int j = 0; j < N; j++)
        {
            if (result[j] != -1)
                res[i][j] = result[j];
            else
            {
                res[i][j] = (i >> ind) & 1;
                ind++;
            }
        }
    }
}

void initialize_vg(NTL::vec_ZZ &vg_, std::map<int, int> c)
{
    static int k = 2;
    vg_.SetLength(r);
    int counter = 0;
    for (int i = -k; i < k + 1; ++i)
    {
        for (int j = 0; j < c[i]; ++j)
        {
            vg_[counter++] = i;
        }
    }

    permuteVector(vg_);

}

NTL::mat_ZZ GSO(const NTL::mat_ZZ& input) {
    NTL::mat_RR mu;
    NTL::vec_RR c;
    ComputeGS(input, mu, c);
    return input;
}

void cyclic_shift_row(NTL::mat_ZZ& mat, long rowIndex, long shiftAmount)
{
    long numCols = mat.NumCols();
    shiftAmount = shiftAmount % numCols; // Ensure the shift amount is within the range of the number of columns

    // Perform the cyclic shift operation
    NTL::vec_ZZ& row = mat[rowIndex];
    NTL::vec_ZZ shiftedRow;
    shiftedRow.SetLength(numCols);

    for (long i = 0; i < numCols; i++)
    {
        long newIndex = (i + shiftAmount) % numCols;
        shiftedRow[newIndex] = row[i];
    }

    // Update the row in the matrix with the shifted values
    row = shiftedRow;
}

NTL::vec_ZZ NP(const NTL::mat_ZZ& B, const NTL::mat_ZZ& B_gs, const NTL::vec_ZZ& t) {

    NTL::vec_ZZ v = t;
    for (int j = N-1; j >= 0; j--) 
    {
        NTL::ZZ B_gs_t_innerProd;
        InnerProduct(B_gs_t_innerProd, B_gs[j], v);

        NTL::ZZ B_gs_self_innerProd;
        InnerProduct(B_gs_self_innerProd, B_gs[j], B_gs[j]);

        NTL::ZZ c =  B_gs_t_innerProd / B_gs_self_innerProd;

        v = v - c * B[j];
    }
    return t-v;
}

NTL::vec_ZZ get_key_vector(NTL::vec_ZZ &vl)
{
    static const int y = 2;

    NTL::vec_ZZ res;

    NTL::mat_ZZ A_z, B_z;
    A_z.SetDims(1, N);
    B_z.SetDims(1, N);

    // set A's boxes
    vectors_for_A(y, vl, A_z);

    int rowsNum = A_z.NumRows();
    res.SetLength(rowsNum << 1);

    for (int q = 0; q < rowsNum; ++q)
    {
        for (int i = 0; i < N; ++i)
        {
            res[q] = res[q] << 1 | A_z[q][i];
        }
    }

    // set B's boxes
    vectors_for_A(y, -vl, B_z);
    rowsNum = B_z.NumRows();
    for (int q = 0; q < rowsNum; ++q)
    {
        for (int i = 0; i < N; ++i)
        {
            res[q + rowsNum] = res[q + rowsNum] << 1 | B_z[q][i];
        }
    }

    std::partial_sort_copy(res.begin(), res.end(), res.begin(), res.end());

    return res;
}

void hybrid_attack(KeyGen & key) {
    NTL::mat_ZZ B, B_reduced, B_gs, A;
    B = Q * NTL::ident_mat_ZZ(N);
    B_reduced = B;

    std::cout << "beta: " << beta << std::endl;
    std::cout << "delta: " << delta << std::endl;
    //std::cout << B << std::endl;

    std::cout << "matrix reduction..." << std::endl;
    NTL::BKZ_FP(B_reduced, delta, beta);
    B_gs = GSO(B_reduced);

    // infinity norm of vl, vg
    const int y = 2, k = 2;

    // ci entries
    std::map<int, int> c;
    int ostatok = r/N * to_int(key.d2 + key.d2_ + key.d1 + key.d1_);
    c.insert({-2, r/N * to_int(key.d2)}); //red
    c.insert({2, r/N * to_int(key.d2)});   // green
    c.insert({-1, r/N * to_int(key.d1)}); //gray
    c.insert({1, r/N * to_int(key.d1)}); //pink
    c.insert({0, N - ostatok}); //white

    while (true)
    {
        // guess vg'
        NTL::vec_ZZ vg_;
        initialize_vg(vg_, c);

        // rotation matrix
        NTL::mat_ZZ A;
        NTL::ZZ_pX aq_px;
        NTL::ZZX aq_x;
        NTL::conv(aq_px, key.aq_invert);
        NTL::conv(aq_x, aq_px);

        A.SetDims(N, N);
        for (int q = 0; q < N; ++q)
        {
            for (int i = 0; i < N; ++i)
            {
                A[q][i] = aq_x[i];
            }
            cyclic_shift_row(A, q, q);
        }
        A = 2 * A;

        // calculate vl' by The Nearest Plane Algorithm
        NTL::vec_ZZ vl_ = NP(B, B_gs, A * vg_);
        NTL::vec_ZZ keyVg_ = get_key_vector(vl_);

        // create an instance of an engine.
        std::random_device rnd_device;
        // Specify the engine and distribution.
        std::mt19937 mersenne_engine {rnd_device()};
        // Generates random integers
        std::uniform_int_distribution<int> dist {-k, k};
        auto gen = [&dist, &mersenne_engine](){
            return dist(mersenne_engine);
        };

        for (int q = 0; q < std::min((long)std::pow(5, N), 10000l); ++q)
        {
            NTL::vec_ZZ vg__;
            vg__.SetLength(r);
            std::generate(std::begin(vg__), std::end(vg__), gen);

            // same vector -> not collision
            if (vg_ == vg__)
            {
                continue;
            }

            NTL::vec_ZZ vl__ = NP(B, B_gs, A * vg__);
            NTL::vec_ZZ keyVg__ = get_key_vector(vl__);

            // collision
            if (keyVg_ == keyVg__)
            {
                NTL::vec_ZZ vg = vg_ + vg__;

                NTL::ZZX vgX;
                NTL::conv(vgX, vg);
                std::cout << vg << std::endl;

                // check found vector
                if (vgX == key.g)
                {

                    std::cout << "Ok" << std::endl;
                    std::cout << "time: " << (float)(clock() - attackTime) / CLOCKS_PER_SEC << " s" << std::endl;
                    std::cout << "vg: " << vg << std::endl;
                    std::cout << "vg': " << vg_ << std::endl;
                    std::cout << "vg'': " << vg__ << std::endl;
                    exit(0);
                }
            }
        }
    }
}



int main(int argc, char **argv)
{

    Setup setup;
    Entropy random;
    Sampler sampler(sigma, alpha_rejection, &random);
    KeyGen key(setup, &random);

    attackTime = clock();
    hybrid_attack(key);
    return 0;
}