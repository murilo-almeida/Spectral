#ifndef _Quadrilateral_headers
#define _Quadrilateral_headers
// ****************************************************************************
class Quadrilateral : public Stdel
// ****************************************************************************
{
 public:
  
  Quadrilateral(int=3,int=4);
  ~Quadrilateral();
  void print(FILE *);
  //void make_Phi(const int m,double Phi[]);
  void make_local_matrices();
  void make_gbnmap(int n0,int n1,int n2,int n3,int ng0,int ng1,int ng2,int ng3,
		  int sign0,int sign1,int sign2,int sgn3,int map[],int sgn[] );

 // Virtual functions
# include "virtual.H"

 private:
  
  void gauss_parameters_default();
 
  //    D ----- C
  //    |       |
  //    |       |
  //    A ----- B
  
  const int aresta[4][2] =
  { {0,1},
    {1,2},
    {3,2},
    {0,3}
  };
  // sinal de rotacao da aresta para obter a normal;
  // +1 = pi/2 anti-horario; -1 = pi/2 horario;
  const int sinal_normal[4] = {-1,-1,1,1};
  
  const int face[1][4] = {{0,1,2,3}};//ABCD// Analisar com cuidado!
  
  int ** ind_mode_;
//  int * emapv; //!< border_map_vertives
//  int * emapi; //!< border_map_indices
};
#endif
