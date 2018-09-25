//
//  Particao.cpp
//  
//
//  Created by Murilo Almeida on 10/26/13.
//
//

#include "spectral.h"

void PARTICAO::Aloca_memoria()
{
	ele.resize(nele);
	bor.resize(nbor);
	gho.resize(ngho);
};
