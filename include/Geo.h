//
//  Geo.h
//  SDG
//
//  Created by Murilo Pereira de Almeida on 10/29/13.
//  Copyright (c) 2013 Murilo Almeida. All rights reserved.
//

#ifndef SDG_Geo_h
#define SDG_Geo_h


/*****************************************************************************/
struct Vertice /*! Vertice = Ponto geometrico */
/*****************************************************************************/
{
  double x,y,z;
  int part_num;
	void set_Vertice(double a,double b,double c)
	{
		x=a;
		y=b;
		z=c;
	};
};
/*! \struct Vertice
 * Vertices with coordinates x,y and z,
 * belonging to partition number part_num
 * \param x = coordinate x
 * \param y = coordinate y
 * \param z = coordinate z
 * \param part_num = number of the partition that holds it
 */

/*****************************************************************************/
struct Elemento { /*! Elemento = Dados para definir o elemento */
/*****************************************************************************/
  int ind;
  int tipo;
  int numv;
  std::vector<int>  vert;
  
  Elemento();
  Elemento(int nv) { numv = nv; vert.resize(nv);/* = new int [nv];*/ }
  ~Elemento();//{ delete [] vert; vert=nullptr;}
};
/*! \struct Elemento
 * Dados de entrada dos elementos
 * \param ind = indice do elemento
 * \param tipo = tipo do elemento
 * \param numv = numero de vertices
 * \param vert[] = conjunto com os indices dos vertices do elemento
 */

struct BCRegion {
  int elnum;
  int eltipo;
  int elface;
};

struct BCStruct {
    int tipo;
    int nregions;
  std::vector<BCRegion>  regions;
    
    BCStruct();
    
    BCStruct(int t, int n) {
        tipo = t;
        nregions = n;
      regions.resize(n);// = new BCRegion [n];
    };
  //~BCStruct(){ delete [] regions; regions=nullptr; }
};

#endif
