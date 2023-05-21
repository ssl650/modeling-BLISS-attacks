#include <boost/program_options.hpp>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include "Setup.h"
#include "KeyGen.h"
#include "Sampler.h"
#include "Sign.h"
#include "Verify.h"
#include "Entropy.h"
#include <ctime>
#include <mpfr.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "ansicolor.h"

#define outof(x,n,tol,s) ((x==n)?Color::FG_XRED:(x>=n-tol?Color::FG_XBLUE:Color::FG_DEFAULT)) << s << Color::FG_DEFAULT

using namespace Eigen;
namespace po = boost::program_options;

typedef Matrix<long, Dynamic, Dynamic> MatrixXl;
typedef Matrix<long, Dynamic, 1>       VectorXl;

long gradasc(VectorXd& v, const MatrixXd& m, const double norm) {
    ArrayXd  a    = ArrayXd::Ones(m.rows());
    VectorXd grad(v.size()), newv(v.size()), dv(v.size());
    
    const double c = 2./(sigma * sigma), phimin = 0.005, nu = 0.8;
    double phi, phi0 = 0.25;
    long iter = 0;

    double l, newl = -INFINITY;

    v = norm * (a.matrix().transpose() * m).transpose().normalized();
    a = (m * v).array();
    a = exp(-c * a);
    l = -log(1 + a).sum();

    for(iter = 0; iter < 20; iter++) {
	a = c / (1 + a.inverse());

	grad = a.matrix().transpose() * m;
	dv   = grad - v.dot(grad) * v / (norm*norm);
	dv   = norm * dv.normalized();

	for(phi = phi0; phi > phimin; phi *= nu) {
	    newv = cos(phi) * v + sin(phi) * dv;
	    a    = (m * newv).array();
	    a    = exp(-c * a);
	    newl = -log(1 + a).sum();

	    if(newl > l + 0.5 * phi * dv.dot(grad))
		break;
	}
	if(newl > l) {
	    v = newv;
	    l = newl;
	}
	if(phi < phi0)
	    phi0 = phi / nu;

	if(phi <= phimin)
	    break;
    }

    return iter;
}


int main(int argc, char *argv[]) {
  long i;

  long nbSign;
  bool compress = false;

  std::string outfile = "stdout";

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "help message")
    ("compress", "generate compressed signatures")
    ("out,o", po::value<std::string>(), "output file")
    ("num,n", po::value<long>(&nbSign)->default_value(1000000), "number of signatures to generate")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return EXIT_SUCCESS;
  }

  if (vm.count("compress"))
    compress = true;

  if (vm.count("out")) {
    outfile.assign(vm["out"].as<std::string>());
    if(std::freopen(outfile.c_str(), "w", stdout) == NULL) {
      std::cerr << "Could not open file '" << outfile << "'."
        << std::endl;
      return EXIT_FAILURE;
    }
  }

  Setup setup;
  Entropy random;
  Sampler sampler(sigma, alpha_rejection, &random);

  std::cerr << "BLISS-" << CLASS << std::endl;

  KeyGen key(setup, &random);
  Sign sign(setup, &sampler, &random);
  Verify verify;

  std::string message = "Hello World!";

  int sign_flip;

  long indc[kappa], z1[N], z2[N], w[2*N], s[2*N];

  Map<VectorXl> vecw1(w,N), vecw2(w+N,N);

  for(i=0; i<N; i++) {
    s[i]   = NTL::to_long(NTL::coeff(key.sk.s1, i));
    s[i+N] = NTL::to_long(NTL::coeff(key.sk.s2, i));
  }

  Map<VectorXl> vecs1(s,N), vecs2(s+N,N), vecz1(z1,N);
  Map<VectorXl> vecs(s,2*N), vecw(w,2*N);

  VectorXd      vecs1d = vecs1.template cast<double>();
  VectorXd      vecs2d = vecs2.template cast<double>();
  VectorXd      vecsd  = vecs.template cast<double>();
  double        s1norm = vecs1d.norm();

  VectorXl	vecw1acc = VectorXl::Zero(N);
  VectorXl	vecw2acc = VectorXl::Zero(N);
  VectorXl      vecv1(N), vecv2(N), vecv(2*N);
  VectorXd      vecv1d(N), vecv2d(N), vecvd(2*N);

  MatrixXd	matw   = MatrixXd::Zero(nbSign, 2*N);

  long		partmatch = -1, fullmatch = -1;

  for(i=0; i<nbSign; i++) {
    sign.signMessage(key.pk, key.sk, message, sign_flip);
    for(int k = 0; k < kappa; k++)
    indc[k] = sign.signOutput.indicesC[k];
    for(int j = 0; j < N; j++) {
    z1[j] = sign.signOutput.z1[j];
    z2[j] = compress ? (sign.signOutput.z2Carry[j]<<dropped_bits) :
	sign.signOutput.z2[j];
    }
    for(int j = 0; j < N; j++) {
	long z1c = 0, z2c = 0;
	for(int k = 0; k < kappa; k++) {
	    long m = j + indc[k];
	    assert( 0 <= m && m < 2*N );
	    if( m < N ) {
	    z1c += z1[m];
	    z2c += z2[m];
	    }
	    else {
	    z1c -= z1[m-N];
	    z2c -= z2[m-N];
	    }
	}
	w[j]  = z1c;
	w[j+N]= z2c;
    }
    
    vecw1acc += sign_flip * vecw1;
    vecw2acc += sign_flip * vecw2;

    matw.row(i) = sign_flip * vecw.template cast<double>();

    if(i%10000 == 0) {
      long matchs1;

      gradasc(vecv1d, matw.topLeftCorner(i+1, N), s1norm);
  
      for(int j = 0; j < N; j++) {
  	vecv1d[j] = std::round(vecv1d[j]);
      }
      vecv1 = vecv1d.template cast<long>();
  
      matchs1 = 0;
  
      for(int j = 0; j < N; j++) {
  	if(vecv1[j] == s[j]) matchs1++;
      }

      std::cerr << outof(matchs1,N,8,"*");
      std::cerr.flush();
      
      if(matchs1 >= N-8 && partmatch < 0)
	  partmatch = i;
      if(matchs1 >= N) {
	  fullmatch = i;
	  break;
      }
    }
  }
  std::cerr << std::endl;

  std::cout << "BLISS-" << CLASS << std::endl;
  std::cout << "Full match: " << fullmatch << std::endl;
  std::cout << "Part match: " << partmatch << std::endl;

  return 0;
}

