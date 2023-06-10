#include "Parameters.h"
#include "Setup.h"
#include "KeyGen.h"
#include "Sign.h"
#include "Verify.h"
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
#define beta 20 // block size for BKZ
#define m (2 * N) // lattice dimension
#define r (N / 2) // meet-in-the-middle dimension
#define delta 0.99
#define prune 10 //???

// Function to permute an NTL vector
void permuteVector(NTL::vec_ZZ& vector)
{
    // Perform permutation
    std::random_device rd;
    std::mt19937 gen(rd());
    shuffle(vector.begin(), vector.end(), gen);
}

void initialize_vg_(NTL::vec_ZZ &vg_, std::map<int, int> c)
{
    static int k = 2;
    vg_.SetLength(r);

    int counter = 0;
    for (int i = -2; i < k + 1; ++i)
    {
        for (int j = 0; j < c[i]; ++j)
        {
            vg_[counter++] = i;
        }
    }

    std::cout << "Default vg':";
    std::cout << vg_ << std::endl;

    permuteVector(vg_);

    std::cout << "Permuted vg':";
    std::cout << vg_ << std::endl;
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

    for (int j = N - 1; j >= 0; j--) {
        NTL::ZZ B_gs_t_innerProd;
        InnerProduct(B_gs_t_innerProd, B_gs[j], v);

        NTL::ZZ B_gs_self_innerProd;
        InnerProduct(B_gs_self_innerProd, B_gs[j], B_gs[j]);

        NTL::ZZ c =  B_gs_t_innerProd / B_gs_self_innerProd;

        v = v - c * B[j];
    }
    return t-v;
}


void hybrid_attack(KeyGen & key) {
    NTL::mat_ZZ B, B_gs, A;
    B = Q * NTL::ident_mat_ZZ(N);
    B_gs = B;

    std::cout << "beta: " << beta << std::endl;
    std::cout << "delta: " << delta << std::endl;
    //std::cout << B << std::endl;

    std::cout << "matrix reduction..." << std::endl;
    NTL::BKZ_FP(B_gs, delta, beta, prune);

    // infinity norm of vl, vg
    const int y = 2, k = 2;

    // ci entries
    std::map<int, int> c;
    c.insert({-2, 29});
    c.insert({2, 29});
    c.insert({-1, 99});
    c.insert({1, 99});
    c.insert({0, 256});

    while (true)
    {
        // guess vg'
        NTL::vec_ZZ vg_;
        initialize_vg_(vg_, c);

        // rotation matrix
        NTL::mat_ZZ A;
        NTL::ZZ_pX aq_px;
        NTL::ZZX aq_x;
        NTL::conv(aq_px, key.aq);
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
        std::cout << "vl':" << std::endl;
        std::cout << vl_ << std::endl;
    }
}



int main(int argc, char **argv)
{

    Setup setup;
    Entropy random;
    KeyGen key(setup, &random);

    std::string message = "IBKS SILA!";

    hybrid_attack(key);
    return 0;
}