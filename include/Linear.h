#ifndef _Linear_headers
#define _Linear_headers
class Linear : public Stdel
{
 public:
  
  Linear(int p0=1,int q0=2);
  ~Linear();
  void print(FILE *);
  //void make_Phi(const int m,double Phi[]);
  void make_local_matrices();
  void make_gbnmap(int n0,int n1,int ng0,int map[],int sgn[] );

	  
  // Virtual functions
 # include "virtual.H"

 private:

  const int aresta[1][2] = {{0,1}};
  
  // sinal para multiplicar (v1 - v0) e obter a normal;
  const int sinal_normal[2] = {-1,1};
  void gauss_parameters_default();
  int * ind_mode_;
  //  int * emapv; //!< border_map_vertives
  //  int * emapi; //!< border_map_indices
};
#endif
