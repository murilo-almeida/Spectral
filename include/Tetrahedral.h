#ifndef _Tetrahedral_headers
#define _Tetrahedral_headers
// class Tetrahedral
// ****************************************************************************
class Tetrahedral : public Stdel 
// ****************************************************************************
{
 public:
  
  Tetrahedral(int p0=3,int q0=4);
  ~Tetrahedral();
  void print(FILE *);
  void make_local_matrices();
 // void ordenar4(int n[], int f_new[]);
  void teste_face(int a,int b,int c, int& sign, int& ng,
		  int& NG,int& NF,int Face[],int Fng[]);
  void map_aresta(int imin,int imax,int sgn0,int ng0,int map[],int sgn[]);
  void map_face(int & imin,int P0,int P1,int sign0,int ng0,int map[],int sgn[]);
  //int show_vtk_type() const {return vtk_type;};
  void localFaceModeMap(const int fnum,const int _P,int trimap[]);
  
  /*
  void Dirichlet(const int face,
                              const Vertice vert[],
                              const int vert_map[],
                              const int gbnmap[],
                              const int sgn[],
                              int bflag[],
                              double Xbc[],
                              double (*func)(double,double,double),
                              const int varn, const int NFields);
*/
  //void printStdel() const;

// virtual functions
# include "virtual.H"
 
 private: 
  int a_CD; // first index of border CD (colapsed border)
  void gauss_parameters_default();
  
  /*       3  <--  colapsado
          /|\
         / | \
        /  |  \
       0 --|---2  <--  colapsado
        \  |  /
         \ | /
          \|/
           1
   */
  const int aresta[6][2] =
  { {0,1}, // AB = 0
    {0,2}, // AC = 1
    {0,3}, // AD = 2
    {1,2}, // BC = 3
    {1,3}, // BD = 4
    {2,3}  // CD = 5 // aresta colapsada : sempre a ultima
  };
  // direcoes nas arestas: dir0
  const int ad0[6] = {0,1,2,1,2,2}; // ad0  = ver[1] - 1;
  // Valores nas outras coordenadas
  // Valores na direcao 1 (dir1): implicito que dir1 = (dir0 + 1) % 3
  const int av1[6] = {0,0,0,0,1,0};
  // Valores na direcao 2 (dir2): implicito que dir2 = (dir0 + 2) % 3
  const int av2[6] = {0,0,0,1,0,1};
 
  const int face_tipo[4] = {2,2,2,2}; //Triangulos
  const int nvf[4] = {3,3,3,3};
  const int face[4][3] = //  sequencia de faces adotada no Gambit
  {
    {0,1,2}, // ABC
    {0,1,3}, // ABD
    {0,2,3}, // ACD
    {1,2,3}  // BCD
  };
  const int face_aresta[4][3] = //  arestas que compoem a face
  {
    {0,3,1}, // ABC = AB, BC, AC
    {0,4,2}, // ABD = AB, BD, AD
    {1,5,2}, // ACD = AC, CD, AD
    {3,5,4}  // BCD = BC, CD, BD
  };

  // direcoes das faces
  const int fd0[4] = {0,0,1,1}; // direcao 0 = ver[1] - 1;
  const int fd1[4] = {1,2,2,2}; // direcao 1 = ver[2] - 1;  dir 2 = 3 - dir0 - dir1
  const int fv2[4] = {0,0,0,1}; // valor da coordenada na direcao normal ( 0 ou 1)
  // sinal de multiplicacao do vetor (l0 X l1) para obter a normal externa;
  const int sinal_normal[4] = {-1,1,-1,1};
  
  int *** ind_mode_;
//  int * emapv; //!< border_map_vertives
//  int * emapi; //!< border_map_indices

};
#endif
