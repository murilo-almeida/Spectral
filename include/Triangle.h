#ifndef _Triangle_headers
#define _Triangle_headers
// ************************************************
class Triangle : public Stdel
{
/* 
     2
    / \
   /   \
  0-----1  
 */
  
 public:
  
  Triangle(int p0=3,int q0=4);
  ~Triangle();
  void print(FILE *);
  //void make_Phi(const int m,double Phi[]);
  void make_local_matrices();
  void make_gbnmap(int n0,int n1,int n2,int ng0,int ng1,int ng2,
                   int sign0,int sign1,int sign2,int map[],int sgn[] );
 
 	  
  // Virtual functions
# include "virtual.H"

 private:
  
 const int aresta[3][2] =
  { {0,1},
    {1,2},
    {0,2}
  };
  // sinal de rotacao da aresta para obter a normal;
  // +1 = pi/2 anti-horario; -1 = pi/2 horario;
  const int sinal_normal[3] = {-1,-1,1};
  
  const int face[1][3] ={{0,1,2}};

  void gauss_parameters_default();
  int ** ind_mode_;
  //  int * emapv; //!< border_map_vertives
  //  int * emapi; //!< border_map_indices
};
#endif
