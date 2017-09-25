#ifndef _ED_Prob_headers
#define _ED_Prob_headers
// ****************************************************************************
// Classe ED : Derivada de GeProb (Generic Problem)
// Encapsula os dados para resolucao do problema de elastodinamica
// ****************************************************************************

//enum var_nomes {ux,uy} ; // Esta ordem influencia na leitura dos parametros

#include "PhElem.hpp"
#include "GeProb.hpp"

// Esta linha eh importante: MyElem eh usado em varios lugares
typedef PhElem<2> MyElem;

// GeProb com elemento MyElemen=PhElem<2> = e duas variaveis
// o segundo 2 abaixo equivale a dois espacos interpolantes

class ED_Prob : public GeProb<MyElem,2,2>
{
 public:

  ED_Prob(Epetra_Comm& comm);
  ~ED_Prob();
  void preamble(char * str);
  void Driver(char * argv = 0);
  void SolveElast2D_dyn(FILE * fout, double lambda, double mu);

  private:

  const int N_VAR = 2;
  const int NFIELDS =2;
}
