//
//  Particao.h
//  
//
//  Created by Murilo Almeida on 10/26/13.
//
//

#ifndef ____Particao__
#define ____Particao__

#include <iostream>
/*****************************************************************************/
struct PARTICAO  /*! PARTICAO = Dados contendo vertices, elementos e bordas
									da particao */
/*****************************************************************************/
{
	PARTICAO(){};
  ~PARTICAO(){};
	void Aloca_memoria();
  int nele,nbor,ngho;
	std::vector<int>  ele;
	std::vector<int>  bor;
	std::vector<int>  gho;
};
/*! \struct PARTICAO
 * \param nele = numero de elementos
 * \param nbor = numero de bordas
 * \param ngho = numero de elementos ghosts
 * \param ele = arrray dos elementos
 * \param bor = array das bordas
 * \param gho = array dos elementos ghosts
 */

#endif /* defined(____Particao__) */
