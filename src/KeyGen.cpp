/*

Copyright or © or Copr. Leo Ducas and Tancrede Lepoint.

Leo.Ducas@ens.fr and Tancrede.Lepoint@ens.fr

This software is a computer program whose purpose is to provide to the
research community a proof-of-concept implementation of the BLISS
digital signature scheme of Ducas, Durmus, Lepoint and Lyubashevsky
appeared at Crypto 2013.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/

#include "KeyGen.h"
#include <NTL/ZZ_pE.h>
#include <iostream>
#include <fstream>

#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <map>

static MYFFT fftKeygen;

#define DEBUG 1

/*
	KeyGen class constructor
*/
KeyGen::KeyGen(const Setup setup, Entropy* random) {

	NTL::ZZ_p::init(NTL::to_ZZ(Q));
	NTL::ZZ_pX phi;
	conv(phi, setup.get_phi());
	NTL::ZZ_pE::init(phi);

	/* Creation private key */
	long i;
	unsigned long randomLong;
	unsigned long randomBit;

	long nbNonZero = setup.get_nb_nonzero();
	long nbNonZero_ = nbNonZero / 2;
    nbNonZero -= nbNonZero_;
	long nbNonZero2 = setup.get_nb_nonzero2();
	long nbNonZero2_ = nbNonZero2 / 2;
    nbNonZero2 -= nbNonZero2_;

    d1 = nbNonZero;
    d1_ = nbNonZero_;
    d2 = nbNonZero2;
    d2_ = nbNonZero2_;

    long normNkS = setup.get_bound_NkS();
#if DEBUG
	std::cout << "nbNonZero = " << setup.get_nb_nonzero() << std::endl;
	std::cout << "nbNonZero2 = " << setup.get_nb_nonzero2() << std::endl;
#endif

	startkg:
	NTL::clear(sk.s1);
	NTL::clear(sk.s2);
	NTL::clear(g);

	i=0;
	while (i < nbNonZero)
	{
		randomLong = random->getRandomLong();
		randomLong = randomLong&(N-1);

		if (NTL::to_long(NTL::coeff(sk.s1, randomLong)) == 0)
		{
			NTL::SetCoeff(sk.s1, randomLong, 1);
			i++;
		}
	}

    i=0;
    while (i < nbNonZero_)
    {
        randomLong = random->getRandomLong();
        randomLong = randomLong&(N-1);

        if (NTL::to_long(NTL::coeff(sk.s1, randomLong)) == 0)
        {
            NTL::SetCoeff(sk.s1, randomLong, -1);
            i++;
        }
    }

	i=0;
	while (i < nbNonZero2)
	{
		randomLong = random->getRandomLong();
		randomLong = randomLong&(N-1);

		if (NTL::to_long(NTL::coeff(sk.s1, randomLong)) == 0)
		{
			NTL::SetCoeff(sk.s1, randomLong, 2);
			i++;
		}
	}

    i=0;
    while (i < nbNonZero2_)
    {
        randomLong = random->getRandomLong();
        randomLong = randomLong&(N-1);

        if (NTL::to_long(NTL::coeff(sk.s1, randomLong)) == 0)
        {
            NTL::SetCoeff(sk.s1, randomLong, -2);
            i++;
        }
    }

	i=0;
	while (i < nbNonZero)
	{
		randomLong = random->getRandomLong();
		randomLong = randomLong&(N-1);

		if (NTL::to_long(NTL::coeff(sk.s2, randomLong)) == 0)
		{
			NTL::SetCoeff(sk.s2, randomLong, 1);
			NTL::SetCoeff(g, randomLong, 1);
			i++;
		}
	}

    i=0;
    while (i < nbNonZero_)
    {
        randomLong = random->getRandomLong();
        randomLong = randomLong&(N-1);

        if (NTL::to_long(NTL::coeff(sk.s2, randomLong)) == 0)
        {
            NTL::SetCoeff(sk.s2, randomLong, -1);
            NTL::SetCoeff(g, randomLong, -1);
            i++;
        }
    }

	i=0;
	while (i < nbNonZero2)
	{
		randomLong = random->getRandomLong();
		randomBit = (randomLong&N)>>log2N;
		randomLong = randomLong&(N-1);

		if (NTL::to_long(NTL::coeff(sk.s2, randomLong)) == 0)
		{
			NTL::SetCoeff(sk.s2, randomLong, 2);
            NTL::SetCoeff(g, randomLong, 2);
			i++;
		}
	}

    i=0;
    while (i < nbNonZero2_)
    {
        randomLong = random->getRandomLong();
        randomBit = (randomLong&N)>>log2N;
        randomLong = randomLong&(N-1);

        if (NTL::to_long(NTL::coeff(sk.s2, randomLong)) == 0)
        {
            NTL::SetCoeff(sk.s2, randomLong, -2);
            NTL::SetCoeff(g, randomLong, -2);
            i++;
        }
    }

#if DEBUG
	std::cout << "f =" << sk.s1 << std::endl;
	std::cout << "g =" << sk.s2 << std::endl;
	std::cout << "g =" << g << std::endl;
#endif

	sk.ls1 = new unsigned char[2*N-1];
	sk.ls2 = new unsigned char[2*N-1];
	sk.offset = N-1;
	for (i=1; i<N; i++)
	{
	  sk.ls1[sk.offset+i] = (unsigned char) (NTL::to_long(coeff(sk.s1, i)) +2);
	  sk.ls1[i-1] = (unsigned char) (-NTL::to_long(coeff(sk.s1, i)) +2);
	  sk.ls2[sk.offset+i] = (unsigned char) NTL::to_long(coeff(sk.s2, i) +2);
	  sk.ls2[i-1] = (unsigned char) (-NTL::to_long(coeff(sk.s2, i)) +2);
	}

	sk.ls1[sk.offset+0] = (unsigned char) (NTL::to_long(coeff(sk.s1, 0)) +2);
	sk.ls2[sk.offset+0] = (unsigned char) (NTL::to_long(coeff(sk.s2, 0)) +2);

    // s2 = 2 * g + 1
	sk.s2 = sk.s2*2+1;
#if DEBUG
	std::cout << "2g+1 =" << sk.s2 << std::endl;

	std::cout << "Norm^2 s1||s2: " << norm2(sk.s1)+norm2(sk.s2) << std::endl;
#endif

	long keyNkS = KeyGen::normNkS();
#if DEBUG
	std::cout << "keyNkS=" << keyNkS << "/" << normNkS << std::endl;
#endif
    // TODO: find out about parameters requirements
	if (keyNkS>normNkS) {
		goto startkg;
	}
	std::cout << "Secret key generated..." << std::endl;
	/* Creation public key */
	NTL::ZZ_pE aq;
	NTL::ZZ_pX pX;
	conv(pX, sk.s1);
	NTL::conv(aq, pX);
	NTL::inv(aq, aq); // aq=1/f mod phi
	NTL::ZZ_pE tmp;
	conv(pX, -(sk.s2));
	NTL::conv(tmp, pX);
	NTL::mul(aq, aq, tmp); // aq_invert = -(2g+1)/f

	//std::cout << "aq_invert=" << aq_invert << std::endl;
	/*pX = rep(aq);

	pk.a1 = new long[2*N-1];
	pk.offset = N-1;
	for (i=1; i<N; i++)
	{
		pk.a1[pk.offset+i] = 2*bmodQ(NTL::to_long(NTL::rep(coeff(pX, i))));
		pk.a1[i-1] = -pk.a1[pk.offset+i];
	}
	pk.a1[pk.offset+0] = 2*bmodQ(NTL::to_long(NTL::rep(coeff(pX, 0))));

	pk.a2 = Q+2;

#if DEBUG
	std::cout << "pk=[";
	for (i=0; i<N; i++)
		std::cout << pk.a1[pk.offset+i] << ", ";
	std::cout << "||" << pk.a2 << "]" << std::endl;
#endif

	std::cout << "Public key generated..." << std::endl;

	pk.a_fft = new long[N];
	for (i=0; i<N; i++) {
		pk.a_fft[i]=bmodQ((pk.a1[pk.offset+i]/2)*W[i]);
	}

	fftKeygen.direct(pk.a_fft);
	for (i=0; i<N; i++)
		pk.a_fft[i] = bmodQ(pk.a_fft[i]);
	pk.modulus = 2*Q;*/

    NTL::clear(this->aq_invert);
    NTL::ZZ_pX tmp_aq_px;
    NTL::ZZ_pE tmp_aq_pe;
    NTL::conv(tmp_aq_px, sk.s2);
    NTL::conv(this->aq_invert, tmp_aq_px);
    NTL::inv(this->aq_invert, this->aq_invert);
    NTL::conv(tmp_aq_px, sk.s1);
    NTL::conv(tmp_aq_pe, tmp_aq_px);
    NTL::mul(this->aq_invert, this->aq_invert, tmp_aq_pe);

#ifdef DEBUG
    std::cout << "aq_invert: " << this->aq_invert << std::endl;
    std::cout << "a^-1 * a == 1: " << ((this->aq_invert * -aq == 1) ? "true" : "false") << std::endl;
#endif
}

/*
	KeyGen class destructor
*/
KeyGen::~KeyGen() {
	delete[] sk.ls1;
	delete[] sk.ls2;
	delete[] pk.a1;
	delete[] pk.a_fft;
}

/*
	(square of) L2-norm of a
*/
unsigned long KeyGen::norm2(NTL::ZZX& a)
{
	unsigned long r=0;
	long c;
	for (long i=0; i<=NTL::deg(a); i++)
	{
		c = NTL::to_long(NTL::coeff(a, i));
		r += c*c;
	}
	return r;
}

/*
	Inner product <l1, l2>
*/
long innerproduct(std::vector<long> l1, std::vector<long> l2)
{
	long result = 0, i;
	for (i=0; i<N; i++)
		result += l1[i]*l2[i];
	return result;
}

/*
	Rotation (multiplication by x^i mod x^n+1)
*/
std::vector<long> rotate(std::vector<long> l, int i)
{
	long j, k=0;
	std::vector<long> result(N);

	for (j=N-i; j<N; j++)
	{
		result[k]=-l[j];
		k++;
	}
	for (j=0; j<N-i; j++)
	{
		result[k]=l[j];
		k++;
	}
	return result;
}

/*
	Norm N_kappa(S)
*/
long KeyGen::normNkS()
{
	long i,j;

	std::vector<long> s1(N), s2(N);
	for (i=0; i<N; i++)
	{
		s1[i] = sk.ls1[sk.offset+i]-2;
		s2[i] = 2*(sk.ls2[sk.offset+i]-2);
	}
	s2[0] = s2[0]+1;
	std::vector<long> gramVec(N);
	std::vector<long> gramVecRot(N);
	std::vector< std::vector<long> > matrixGram(N, std::vector<long>(N));

	for (i=0; i<N; i++)
	{
		gramVec[i] = innerproduct(s1, rotate(s1,i))+innerproduct(s2, rotate(s2,i));
	}

	for (i=0; i<N; i++)
	{
		gramVecRot = rotate(gramVec,i);
		for (j=0; j<N; j++)
			matrixGram[i][j]=gramVecRot[j];
	}

	for (i=0; i<N; i++)
		std::sort (matrixGram[i].begin(), matrixGram[i].begin()+N);

	std::vector<long> maxis(N);
	for (i=0; i<N; i++)
	{
		maxis[i] = 0;
		for (j=0; j<kappa; j++)
		{
			maxis[i] += matrixGram[i][N-1-j];
		}
	}
	std::sort(maxis.begin(), maxis.begin()+N);
	long result=0;
	for (i=0; i<kappa; i++)
		result += maxis[N-1-i];
	return result;
}
