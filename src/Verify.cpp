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

#include "Verify.h"

#include <fstream>
#include <openssl/sha.h>
#include <iomanip>
#include <cstring>

MYFFT fftVerify;

/*
	Drop dropped_bits bits in the vector v, and store it in result

	WARNING: should be same code as in Sign.cpp
*/
void roundDV(long* result, long* v)
{
	long i;
	for (i=0; i<N; i++)
		result[i] = (2*v[i]+(1<<dropped_bits))/(1<<(dropped_bits+1));
}

/*
	Verify the signature
*/
bool Verify::verifyMessage(const struct pubkey& pk, const struct signature& s, const std::string message)
{
	long* z2d = new long[N];
	for (long i=0; i<N; i++)
	{
		z2d[i] = (1<<dropped_bits)*s.z2Carry[i];
	}

	// Verification infinite norm
	if (normInfinite(s.z1)>=BINFTY or normInfinite(z2d)>=BINFTY){
		std::cout << "Problem Infinite Norm" << std::endl;
		return 0;
	}

	// Verification l2-norm
	if (norm2(s.z1)+norm2(z2d) > B_VERIFY){
		std::cout << "Problem L2 Norm" << std::endl;
		return 0;
	}

	unsigned char hash[SHA512_DIGEST_LENGTH];
	char arrayMessage[message.size()+1];
	strncpy(arrayMessage, message.c_str(), message.size());
	sha512(hash, arrayMessage, message.size());

	// New Hash for generateC for random oracle
	char newHash[SHA512_DIGEST_LENGTH+N*2+1];
	for (long i=0; i<SHA512_DIGEST_LENGTH; i++)
		newHash[i]=hash[i];

	// construct a/(q+2)*z_1+q/(q+2)*c+z_2
	long* az = new long[N];
	long i;
	for (i=0; i<N; i++)
	{
		az[i] = bmodQ(s.z1[i]*W[i]);
	}
	fftVerify.direct(az);

	for (i=0; i<N; i++)
	{
		az[i] = bmodQ(az[i]*pk.a_fft[i]);
	}
	fftVerify.inverse(az);

	for (i=0; i<N; i++)
	{
		if (az[i]>Q) az[i] -= Q; // 0<= az[i]<2Q: reduction mod QQQ
		az[i] = bmod2Q(2*az[i]*oneQTwo);
	}

	for (i=0; i<kappa; i++) {
		az[s.indicesC[i]] = bmod2Q(az[s.indicesC[i]]+Q*oneQTwo);
	}

	// Verify if c = H(az+qc, mu)
	long* indices = new long[kappa];

	roundDV(az,az);
	for (i=0; i<N; i++)
	{
		az[i] = (az[i]+s.z2Carry[i]);
		if (az[i]<0)
			az[i] += modp;
		if (az[i]>=modp)
		{
			az[i]-=modp;
		}
	}
	generateC(indices, az, newHash);

	delete[] az;

	for (i=0; i<kappa; i++)
		if (indices[i]!=s.indicesC[i])
			{
				delete[] indices;
				return 0;
			}

	delete[] indices;
	return 1;
}

/*
	Compute infinite norm
*/
long normInfinite (long* z)
{
	long max=0, i;
	for (i=0; i<N; i++)
	{
		if (z[i]>0 && z[i]>max) max = z[i];
		else if (z[i]<0 && z[i]<-max) max = -z[i];
	}
	return max;
}

/*
	Generate Challenge

	WARNING: should be same code as in Sign.cpp
*/
void Verify::generateC(long* indices, long* ay, char* newHash)
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
