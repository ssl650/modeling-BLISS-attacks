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

#ifndef _FFT_H
#define _FFT_H

#include "Parameters.h"

typedef unsigned long uint128_t __attribute__((mode(TI)));

class MYFFT {

private:

	inline void butterfly(long* x, long in1, long in2, long inW);
	void rearrange(long* x);

public:
	MYFFT() {};
	void direct(long* x);
	void inverse(long* x);
	~MYFFT() {};

};

template<typename V>
inline V bmodQ(const V&a)
{
	V c = a;
	if (a<0)
		c += Q;
	c = c-Q*((oneoverQ*c)>>32);
	while (c>=Q){ c-= Q;}
	return c;
}

template<typename V>
inline V bmod2Q(const V&a)
{
	V c = a;
	if (a<0)
		c += 2*Q;
	c = c-2*Q*((oneover2Q*c)>>32);
	while (c>=2*Q){ c-= 2*Q;}
	return c;
}

#endif
