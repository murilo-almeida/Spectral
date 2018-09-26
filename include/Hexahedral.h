#ifndef _Hexahedral_headers
#define _Hexahedral_headers

#include "PhElem.hpp"
// ****************************************************************************
class Hexahedral : public Stdel
// ****************************************************************************
{
 public:
  
  Hexahedral(int=3,int=4);
  ~Hexahedral();
  void print(FILE *);
  //void make_Phi(const int m,double Phi[]);
  void make_local_matrices();
  void make_gbnmap(int n0,int n1,int n2,int n3,int ng0,int ng1,int ng2,int ng3,
		  int sign0,int sign1,int sign2,int sgn3,int map[],int sgn[] );
  void superficie_externa(const int Vert_map[],const int & num_local,double & area,double normal[3]);
  void localFaceModeMap(const int fnum, int lmap[]);

 // Virtual functions
# include "virtual.H"

 private:

 void gauss_parameters_default();
 // Ordem local

 //     C ---- D
 //   / |     /|
 // G ------H  |
 // |   A --|- B
 // | /     | /
 // E-------F
 const int aresta [12][2] = 
   { {0,1},
     {0,2},
     {0,4},
     {1,3},
     {1,5},
     {2,3},
     {2,6},
     {3,7},
     {4,5},
     {4,6},
     {5,7},
     {6,7}
   };
 const int face_tipo[6] = {3,3,3,3,3,3}; // Quadrilateros
 const int nvf[6] = {4,4,4,4,4,4};
 const int face[6][4] = 
   { {0,1,3,2},
     {0,1,5,4},
     {0,2,6,4},
     {1,3,7,5},
     {2,3,7,6},
     {4,5,7,6} 
   };
  // direcoes das faces
  const int fd0[6] = {0,0,1,1,0,0}; // direcao 0 = (ver[1] - ver[0])/2;
  const int fd1[6] = {1,2,2,2,2,1}; // direcao 1 = (ver[3] - ver[0])/2;  dir 2 = 3 - dir0 - dir1
  const int fv2[6] = {0,0,0,1,1,1}; // valor da coordenada na direcao normal ( 0 ou 1)
 // incremento no indice dos pontos de Gauss
    int inc[3];
  // sinal de multiplicacao do vetor (l0 X l1) para obter a normal externa;
 const int sinal_normal[6] = {-1,1,-1,1,-1,1};
 int *** ind_mode_;
//  int * emapv; //!< border_map_vertives
//  int * emapi; //!< border_map_indices
};

#endif
