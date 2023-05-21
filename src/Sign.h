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

#ifndef __SIGN_H
#define __SIGN_H

#include "Parameters.h"
#include "Fft.h"
#include "Setup.h"
#include "KeyGen.h"
#include "Sampler.h"

#include<string>

struct signature {
	long* z1;
	long* z2;
	long* z1High;
	long* z1Low;
	long* z2Carry;
	long* indicesC;
};

class Sign {
	private:
		Sampler* sampler;
		Entropy* random;
		long bound_NkS;
	
	public:
		struct signature signOutput;
		Sign(Setup setup, Sampler* s, Entropy* randomGen);
		long signMessage(const struct pubkey& pk, const struct seckey& sk, const std::string message, int& sign_flip);
		void mult_by_c(long* result, unsigned char* ls, bool isG, long offset, long* indices);
		void mult_by_a2_fft(long* ay, long* a_fft, long* y2);
		void computeCarries(long* result, long* ay, long* z2);
		void generateC(long* indices, long* ay, char* hash);
		~Sign();
};

void sha512(unsigned char* hash, char* message, long size);
unsigned long norm2(long* a);
long scalProd(long* z, long* a);

#endif
