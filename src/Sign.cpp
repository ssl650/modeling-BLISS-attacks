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

#include "Sign.h"
#include <fstream>
#include <openssl/sha.h>
#include <iomanip>
#include <cstring>

static MYFFT fftSign;

/*
	Sign class construction
*/
Sign::Sign(Setup setup, Sampler* s, Entropy* randomGen) {
	sampler = s;
	random = randomGen;
	bound_NkS = setup.get_bound_NkS();
	signOutput.z1 = new long[N];
	signOutput.z2 = new long[N];
	signOutput.indicesC = new long[kappa];

	signOutput.z1High = new long[N];
	signOutput.z1Low = new long[N];
	signOutput.z2Carry = new long[N];
}

/*
	Sign class destruction
*/
Sign::~Sign() {
	delete[] signOutput.z1;
	delete[] signOutput.z2;
	delete[] signOutput.indicesC;

	delete[] signOutput.z1High;
	delete[] signOutput.z1Low;
	delete[] signOutput.z2Carry;
}

/*
	Drop dropped_bits bits in the vector v, and store it in result

	WARNING: should be same code as in Verify.cpp
*/
void roundD(long* result, long* v)
{
	long i;
	for (i=0; i<N; i++)
	{
		result[i] = (2*v[i]+(1<<dropped_bits))/(1<<(dropped_bits+1));
		result[i] = result[i]%modp;
	}
}

/*
	Sign a message, and return the number of rejections
*/
long Sign::signMessage(const struct pubkey& pk, const struct seckey& sk, const std::string message, int& sign_flip)
{
	// Define temporary variables
	long i;
	long norm;
	long nb_tests = 0;

	long* sc1;
	long* sc2;
	long* ay;
	long* tmp;
	sc1 = new long[N];
	sc2 = new long[N];
	ay 	= new long[N];
	tmp = new long[N];
	long* z2d = new long[N];

	// Hash the message: hash[] = sha512(message)
	unsigned char hash[SHA512_DIGEST_LENGTH];
	char arrayMessage[message.size()+1];
	strncpy(arrayMessage, message.c_str(), message.size());
	sha512(hash, arrayMessage, message.size());

	// New Hash for generateC for random oracle
	char newHash[SHA512_DIGEST_LENGTH+N*2+1];
	for (i=0; i<SHA512_DIGEST_LENGTH; i++)
		newHash[i]=hash[i];

	// Start Sign
	startSign:
	nb_tests++; // Update nb of tests

	// Gaussian Sampling
	for (i=0; i<N; i++) {
		signOutput.z1[i] = sampler->SamplerGaussian();
		ay[i] = signOutput.z1[i]*W[i]; // Formatting for FFT
		signOutput.z2[i] = sampler->SamplerGaussian();
	}

	// Compute a/(q+2)*y_1 + y_2
	mult_by_a2_fft(ay, pk.a_fft, signOutput.z2);

	// Generate c
	roundD(tmp, ay);
	generateC(signOutput.indicesC, tmp, newHash);

	// Compute s*c
	mult_by_c(sc1, sk.ls1, false, sk.offset, signOutput.indicesC);
	mult_by_c(sc2, sk.ls2, true, sk.offset, signOutput.indicesC);

	// Reject with proba 1/(M*exp(-norm(sc)^2/2*sigma*sigma))
	norm = norm2(sc1)+norm2(sc2);
	if (!sampler->SamplerBerExpM(norm)) goto startSign;

	// Compute z
	if (random->getRandomBit()) {
		for (i=0; i<N; i++) {
			signOutput.z1[i]-=sc1[i];
			signOutput.z2[i]-=sc2[i];
			sign_flip = -1;
		}
	} else {
		for (i=0; i<N; i++) {
			signOutput.z1[i]+=sc1[i];
			signOutput.z2[i]+=sc2[i];
			sign_flip = 1;
		}
	}

	// Reject with probability 1/cosh(<z, sc>/sigma*sigma)
	if (!sampler->SamplerBerCosh(
		scalProd(signOutput.z1, sc1)+scalProd(signOutput.z2, sc2)))
		goto startSign;


	computeCarries(signOutput.z2Carry, ay, signOutput.z2);

	for (long i=0; i<N; i++)
	{
		z2d[i] = signOutput.z2Carry[i] << dropped_bits;
	}

	for (i=0; i<N; i++)
	{
		if (signOutput.z1[i]>=BINFTY || z2d[i]>=BINFTY || signOutput.z1[i]<=-BINFTY || z2d[i]<=-BINFTY) goto startSign;
	}

	if (norm2(signOutput.z1)+norm2(z2d) >= B_VERIFY )
		goto startSign;

	delete[] sc1;
	delete[] sc2;
	delete[] ay;
	delete[] tmp;
	delete[] z2d;

	return nb_tests;
}


/*
	Sha-512 of message into hash
*/
void sha512(unsigned char* hash, char* message, long size)
{
    SHA512_CTX sha512;
    SHA512_Init(&sha512);
    SHA512_Update(&sha512, message, size);
    SHA512_Final(hash, &sha512);
}


/*
	Multiplication of S by c
*/
// Faster ! Multiplication by c without reduction mod 2Q
// We use a unsigned long sum as 8 unsigned char sums
// Valid as long as kappa <= 63
// Any change to the format of the key may break this implementation
void Sign::mult_by_c(long* result, unsigned char* ls, bool isG, long offset, long* indices)
{
  long i, j, corr = -2*kappa;
  unsigned long value;
  unsigned char* char_value  = (unsigned char*) ((void*) (&value));
  unsigned long* tmp;

  for (j=0; j<N; j+=8){
    value= 0;
    for (i=0; i<kappa; i++){
      tmp = (unsigned long*) (ls+(offset+j)-indices[i]);
      value += *tmp;
    }
    result[j+0] = ((long) char_value[0])+corr;
    result[j+1] = ((long) char_value[1])+corr;
    result[j+2] = ((long) char_value[2])+corr;
    result[j+3] = ((long) char_value[3])+corr;
    result[j+4] = ((long) char_value[4])+corr;
    result[j+5] = ((long) char_value[5])+corr;
    result[j+6] = ((long) char_value[6])+corr;
    result[j+7] = ((long) char_value[7])+corr;
  }
  if (isG)
  {
	for (j=0; j<N; j++)
	{
		result[j] = 2*result[j];
	}

	for (i=0; i<kappa; i++)
		result[indices[i]] = result[indices[i]]+1;
  }
}

/*
	Multiplication by a/(q+2) with FFT
*/
void Sign::mult_by_a2_fft(long* ay, long* a_fft, long* y2)
{
	long i;

	fftSign.direct(ay);

	for (i=0; i<N; i++)
		ay[i] = bmodQ(ay[i]*a_fft[i]);

	fftSign.inverse(ay);

	for (i=0; i<N; i++) {
		if (ay[i]>Q) ay[i] -= Q; // 0<= ay[i]<2Q: reduction mod QQQ
		ay[i] = bmod2Q((2*ay[i]*oneQTwo)+y2[i]);
	}
}

/*
	Compute carries of z2
*/
void Sign::computeCarries(long* result, long* ay, long* z2)
{
	long i;

	roundD(result, ay);

	long* tmpC = new long[N];
	for (i=0; i<N; i++)
	{
		tmpC[i] = bmod2Q(ay[i]-z2[i]);
	}
	roundD(tmpC, tmpC);

	for (i=0; i<N; i++)
	{
		long value = result[i]-tmpC[i];
		if (value<-modp/2)
			value += modp;
		if (value>modp/2)
			value -= modp;
		result[i] = value;
	}
	delete[] tmpC;
}

/*
	Generate Challenge

	WARNING: should be same code as in Verify.cpp
*/
void Sign::generateC(long* indices, long* ay, char* newHash)
{
	long i, j;

	for (i=0; i<N; i++)
	{
		for (j=0; j<2; j++)
		{
			newHash[SHA512_DIGEST_LENGTH+i*2+j] = ay[i]&((unsigned char)-1);
			ay[i] >>= 8;
		}
	}

    unsigned char hash[SHA512_DIGEST_LENGTH];
    unsigned char nbRepetitions = 0;
    unsigned char arrayValues[N];

    randomOracle:
    newHash[SHA512_DIGEST_LENGTH+N*2]=nbRepetitions;
    nbRepetitions++;
    sha512(hash, newHash, SHA512_DIGEST_LENGTH+N*2+1);

    for (i=0; i<N; i++)
    	arrayValues[i] = 0;

    unsigned long index;
    j = 0;
 #if CLASS==0 // 8 bits = 1 char

    for (i=0; i<kappa;) {
      index = (long) hash[j];
      if (!arrayValues[index]){
    	indices[i] = index;
    	arrayValues[index]++;
	i++;
      }
      j++;
      if (j>=64) goto randomOracle;
    }

#else // 9 bits
    unsigned long extra_bits;

    extra_bits = *((unsigned long *) (&hash[56]));

    for (i=0; i<kappa;) {
      index = 2*((long) hash[j]) + (extra_bits %2);
      extra_bits = extra_bits >> 1;
      if (!arrayValues[index]){
    	indices[i] = index;
    	arrayValues[index]++;
	i++;
      }
      j++;
      if (j>=56) goto randomOracle;
    }

    #endif

}

/*
	(square of the) L2-Norm of the vector a
*/
unsigned long norm2(long* a)
{
	unsigned long r=0;
	for (long i=0; i<N; i++)
		r += a[i]*a[i];
	return r;
}

/*
	Scalar Product <z,a>
*/
long scalProd(long* z, long* a)
{
	long r=0;
	for (long i=0; i<N; i++)
		r += z[i]*a[i];
	return r;
}
