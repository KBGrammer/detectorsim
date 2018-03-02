#ifndef H_RNG
#define H_RNG


// Random number generator class
class Rng {
 private:
  unsigned long long u,v,w;
  double exc_div, inc_div;

 public:
  Rng(unsigned long long j);
  ~Rng();

  unsigned long long rint64();
  double rdouble();
  float rfloat();
  double rdouble_exc();
  unsigned int rint32();
};

// Gaussian random number generator
class GRng: public Rng {
 private:
  double mu, sigma;

 public:
  GRng(unsigned long long j);
  GRng(double m, double s, unsigned long long j);
  ~GRng();

  double rgauss();
  void set_gauss(double m, double s);
};


#endif
