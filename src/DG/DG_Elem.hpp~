//
//  DG_Elem.hpp
//  SDG
//
//  Created by Murilo Almeida on 30/10/16.
//  Copyright © 2016 Murilo Almeida. All rights reserved.
//

#ifndef DG_Elem_hpp
#define DG_Elem_hpp

#include <stdio.h>
# include "spectral.h"

# include "PhElem.hpp"

// *************************************************************************
class DG_Elem : public PhElem<2>
{
public:
  
  DG_Elem(){ rho=1.0;vetores_iniciados=1;}// 1 (=FALSE)}; //(const int n = 1);
  ~DG_Elem(){/*cout << "destruir DG_Elem\n";*/ };
  
  void inicia_vetores();
  void set_fontes(double sw,double sn){qw=sw; qn=sn;};
  double show_rho(){return rho;};
  double show_Volume(){return Volume;};
  void inicia_funcoes_na_borda( EDGE *  border);
  void finaliza_vetores();
  void inicia_tracos(EDGE *  border);
  void echo_traco(FILE * f_eco = nullptr);

  void Atualizar_valores(FILE * fout = NULL);// NULL eh o valor default
  void set_permeabilidade(const double,const double, const double);
  void set_mass_density(double);
  void set_porosidade(double);
    
  // Funcoes especificas para DG_Problem
  void VolumeIntegrals_UMFPACK(const double Dt,Fluids fls,
                               int & count,
                               int * Ti,int * Tj, double * Tx,
                               double * B,
                               double * = NULL,
                               double * = NULL);
  void VolumeIntegrals(const double Dt,Fluids fls,
                       Teuchos::RCP<Epetra_FECrsMatrix>   A,
                       Teuchos::RCP<Epetra_FEVector> RHS,
                       double * = NULL,
                       double * = NULL);
  void VolumeTracos(const double Dt,Fluids fls,
                    double * = NULL,
                    double * = NULL);
  void VolumeIntegralsT(const double Dt_new,const double Dt_old,int & count,
                        int * Ti,int * Tj, double * Tx,
                        double * B);
  
  void calcula_tracos(Fluids fls);
  // Initial Guess Volume Integrals
  void VolumeIntegrals_IG(Fluids fls, int & count,
                          int * Ti,int * Tj, double * Tx,
                          double * B);
  
  void perm_val(double K[3]);
  
  void Traco_sn(const int & h,double *saida);
  
  void Traco_pw(const int & h,double *saida);
  
  void Traco_Kgrad_pw(const int & h,double ** saida);
  
  void Traco_Kgrad_pc(const int & h,double ** saida);
  
  void Traco_Kgrad_sn(const int & h,double ** saida);
  
  void Traco_phi(const int & lado,const int & ivar,
                 const int & ind,double * saida);
  
  // ***********************************************************************
  void Traco_phi_1(const int & lado,const int & ivar,const int & pos,
                   double * saida)
  {
    int qmax=ptr_stdel[0]->qborder_val();
    for(int q=0;q<qmax;q++)
      saida[q]=TrPhi[ivar][lado][pos][q];
  };
  // ***********************************************************************

  void Traco_grad_phi(const int & lado,const int & ivar,
                      const int & ind,double ** saida);
  void Traco_Kgrad_phi_n(const int & lado,const int & ivar,
                         const int & ind,double * saida);
  
  // void Row(int * NumNz, int ** MapRow); // Nao esta em uso
  // ***********************************************
  // Versao paralela. Dados em vetor
  // ***********************************************
  // void vetorizar_dados(int elm_num, int &count1, int &count2, int *controle, int *conteudo);
  
  //void construir_PhiArray();
  //const double * eval_Phi(const int i);

  
private:
  
  // Variaveis especificas para o caso DG_Problem
  double * sna, * pwa; //!< ponteiros para os valores nos pontos de Gauss da saturacao e pressao atuais
  double rho; //!< densidade de massa ou porosidade
  double porosidade;
  double perm[3]; //!< permeabilidades nas 3 direcoes
  double qn,qw;
  double Volume; //!< Volume do elemento ( = area do elemento para elementos bidimensionais e comprimento para elem. 1D
 
  double     * Jb; //!< Jacobiano nos pontos de Gauss sobre as bordas; usado na integração sobre bordas
  double    ** Mass_sn;
  double    ** Trsn;
  double    ** Trpw;
  double   *** LaplacianoPhi;
  double   *** TrKgrad_sn;
  double   *** TrKgrad_pw;
  double   *** TrKgrad_pc;
  double  **** TrPhi;
  double  **** GradPhi;
  double  **** TrKgradPhi_n;
  double ***** TrGradPhi;
  double     * PhiArray;
  
  //double ** b;
  
};
/*! \class DG_Elem
 */


// ****************************************************************************
// Class DG_Elem ---  Para escoamentos de 2 fluidos imisciveis DG
// ****************************************************************************

#endif /* DG_Elem_hpp */
