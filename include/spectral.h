/** Main header file for the spectral family
 * Using newmat functions for small matrixes and vectors manipulations
 * Using umfpack for sparse matrices
 * Defining the size parameters for the global matrices
*/
#ifndef SPECTRAL_H
#define SPECTRAL_H

# include "MyOptions.h"

// Epetra definitions
#include "MyTrilinos.h"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <stdio.h>
#include <algorithm>
#include <vector>
#include <cstddef>
#include "boost/multi_array.hpp"

#ifdef _NEWMAT
// NEWMAT
# include "include.h"
# include "newmat.h"
#include "newmatio.h"
#endif

extern "C" {
#ifdef _MAC_OS
#include <umfpack.h>
#else
#include <suitesparse/umfpack.h>
#endif
}

using namespace std;

//# include "Mat2.h"
//# include "Matn.h"

struct Field_struct {
  int ordem,P,Q;
};

// tipo numero identico ao gmsh
const int ElemType[6][3]={
  {1,0,0}, //!< 0=ponto
  {2,1,0}, //!< 1=linha
  {3,3,1}, //!< 2=triangulo
  {4,4,1}, //!< 3=quadrilatero
  {4,6,4}, //!< 4=tetraedro
  {8,12,6} //!< 5=hexaedro
};
 /*! \var ElemType
  * Valores usados pelo gmsh.
  * Variavel global contendo: {numv, ne, nf}
  * \param ElemType[][0] = numv = numero de vertices
  * \param ElemType[][1] = ne = numero de edges
  * \param ElemType[][2] = nf = numero de faces
  * \param ElemType[0][] Ponto  tipo=0;
  * \param ElemType[1][] Linear tipo=1;
  * \param ElemType[2][] Triangular tipo=2;
  * \param ElemType[3][] Quadrangular tipo=3;
  * \param ElemType[4][] Tetraedrico tipo=4;
  * \param ElemType[5][] Hexaedrico tipo=5; */

# include "Geo.h" 
/*****************************************************************************/
class Mode  /*! mode  */
/*****************************************************************************/
{
 public:

  Mode(int ap=0,int aq=0,int ar=0) {set_mode(ap,aq,ar);};
  void set_mode(int ap,int aq,int ar){p=ap; q=aq; r=ar;};
  int p_val(){return p;};
  int q_val(){return q;};
  int r_val(){return r;};
  void print(FILE *fout = nullptr){
    if(fout == nullptr) printf("%4d %4d %4d\n",p_val(),q_val(),r_val());
    else
      fprintf(fout,"%4d %4d %4d\n",p_val(),q_val(),r_val());
  };

 private:

  int p,q,r;
};
/*! \class Mode
 * integers specifying the order of the modes in the local x,y and z directions,
 * respectively.
 * \param p = order of the mode in the x-direction
 * \param q = order of the mode in the y-direction
 * \param r = order of the mode in the z-direction
 */
/*****************************************************************************/
struct ARESTA  /*! ARESTA = Aresta do elemento */
/*****************************************************************************/
{
  ARESTA()
  {
    // pdir = std::nullptr_t;
    // sdir = std::nullptr_t;
  };
  ~ARESTA()
  {
    //if(pdir != std::nullptr_t) {
    //delete [] pdir; pdir = nullptr;
    //}
    //if(pdir != std::nullptr_t) {
    //delete [] sdir; sdir = nullptr;
    //}
  };
  int Na,Nb;
};
/*****************************************************************************/
struct FACET  /*! FACET = Aresta em 1D e Face em 2D do elemento */
/*****************************************************************************/
{
  FACET()
  {
    // pdir = std::nullptr_t;
    // sdir = std::nullptr_t;
  };
  ~FACET()
  {
    //if(pdir != std::nullptr_t) {
    //delete [] pdir; pdir = nullptr;
    //}
    //if(pdir != std::nullptr_t) {
    //delete [] sdir; sdir = nullptr;
    //}
  };
  int tipo,num_elem,elemento[2],num_local[2];
};

/*****************************************************************************/
struct EDGE  /*! EDGE = Borda do elemento */
/*****************************************************************************/
{
  EDGE()
  {
    // pdir = std::nullptr_t;
    // sdir = std::nullptr_t;
  };
  ~EDGE()
  {
    //if(pdir != std::nullptr_t) {
    //delete [] pdir; pdir = nullptr;
    //}
    //if(pdir != std::nullptr_t) {
    //delete [] sdir; sdir = nullptr;
    //}
  };
  int tipo,Na,Nb,num_elem,elemento[2],num_local[2],sinal[2];
  double normal[3],comprimento;
  double * pdir;
  double * sdir;
  //  double flux_w, flux_n;
  int gbtrbind[2];
  int part_num;
};

/*! \struct EDGE
 * Arestas das faces, compostas por dois vertices;
 * \param tipo  =-1 inflow;
 *              = 0 no-flow;
 *              =+1 outflow;
 *              = 2 interior;
 * \param Na    = number of the initial mode
 * \param Nb    = number of the final mode
 * \param elemento[2] = numbers of the elements sharing the border
 * \param num_local[2] = side numbers the border assumes in the elements
 * \param sinal[2] = sign the border has in the elements; -1 if the border is
 followed backwards;
 * \param normal[2]= components of the normal vector;
 * \param comprimento = length of the border;
 * \param *pdir = pointer to the array containing the Dirichlet values of pressure
 * \param *sdir = pointer to the array containing the Dirichlet values of saturation
 * \param gbtrbind = global trace border index for the traces on the sides
 * \param part_num = partition number
 * \param flux_w = flux of w through the border
 * \param flux_n = flux of non-wetting fluid through the border
 */
/*****************************************************************************/
struct FACE  /*! FACE = Borda do elemento tridimensional*/
/*****************************************************************************/
{
  FACE()
  {
    // pdir = std::nullptr_t;
    // sdir = std::nullptr_t;
  };
  ~FACE()
  {
    //if(pdir != std::nullptr_t) {
    //delete [] pdir; pdir = nullptr;
    //}
    //if(pdir != std::nullptr_t) {
    //delete [] sdir; sdir = nullptr;
    //}
  };
  int _tipo,_nv,_elemento[2],_face[2];
  std::vector<int> _vertice;
  double _normal[3],_area;
  double * pdir;
  double * sdir;
  //  double flux_w, flux_n;
  //int gbtrbind[2];
  int part_num;
};

/*! \struct FACE
 * Faces dos elementos tridimensionais, compostas por _nv vertices;
 * \param _tipo =;
 *              = ;
 * \param _nv   = number of vertices
 * \param _vertice = vertice number;
 * \param _elemento[2] = numbers of the elements sharing the border
 * \param _face[2] = local face number which the face assumes in the elements
 * \param _normal[3]= components of the normal vector;
 * \param _area = face area;
 * \param *pdir = pointer to the array containing the Dirichlet values of pressure
 * \param *sdir = pointer to the array containing the Dirichlet values of saturation
 * \param part_num = partition number
 * \param flux_w = flux of w through the border
 * \param flux_n = flux of non-wetting fluid through the border
 */

# include "Funcoes_c.h"
# include "Tstruct.h"
# include "Particao.h"
# include "Stdel.h"
# include "Linear.h"
# include "Triangle.h"
# include "Quadrilateral.h"
# include "Tetrahedral.h"
# include "Hexahedral.h"
#endif
