//
//  DG_Elem.cpp
//  SDG
//
//  Created by Murilo Almeida on 30/10/16.
//  Copyright © 2016 Murilo Almeida. All rights reserved.
//

#include "spectral.h"
#include "DG_Elem.h"

// ***************************************

// ****************************************************************************
void DG_Elem::Atualizar_valores(FILE * fileout)
{
  
  //cout << "sat " << endl;
  ptr_stdel[sat]->evalGQ(sna,u0[sat]);
  //cout << "pres " << endl;
  ptr_stdel[pres]->evalGQ(pwa,u0[pres]);
  
  if(fileout != NULL){
    ptr_stdel[sat]->printGQtofile(fileout,sna,pwa,ptvert,Vert_map);
    fflush(fileout);
  }
};

// ****************************************************************************
void DG_Elem::set_permeabilidade(const double x,const double y, const double z)
{
  perm[0]=x;
  perm[1]=y;
  perm[2]=z;
};

// ****************************************************************************
void DG_Elem::set_mass_density(double value)
{
  rho=value;
};

// ***************************************************************************
void DG_Elem::inicia_funcoes_na_borda(EDGE * border)
{
  // Esta procedure permite que se calcule os tracos
  // diretamente das funcoes sem fazer a operacao de
  // fatiamento. Pode ser usada com pontos de Gauss
  // totalmente internos.
  // ***************************************
  // Calcular Phi e seu gradiente na borda *
  // Nao deve usar eval_Phi() pois essa    *
  // usa os pontos de Gauss do elemento,   *
  // que podem nao existir nas bordas      *
  // ***************************************
  int h,i,j,pos;
  double n_e[ndim];
  
  for(i=0;i<NumLocalVars;i++) {//i= variavel
    
    int nn     = numn[i];
    int qmax   = ptr_stdel[i]->qborder_val();
    int NGQP   = ptr_stdel[i]->NGQP_val();
    double phi[NGQP];
    
    //	double TP[nborder][nn][qmax];
    
    double *** TP = new double ** [numborders];
    for(h=0;h<numborders; h++) {
      TP[h] = new double * [nn];
      for(j=0;j<nn;j++){
        TP[h][j] = new double [qmax];
      }
    }
    
    //int s0 = nn*ndim*qmax;
    //int s1 =    ndim*qmax;
    
    // No elemento
    for(j=0;j<nn;j++){
      ptr_stdel[i]->eval_Phi(j,phi);  // <-- No elemento
      // Calcular as derivadas de Phi_j
      ptr_stdel[i]->Gradiente(GradPhi[i][j],phi,ptvert,Vert_map);
    }
    
    // *************************************************
    // Calcula os Laplacianos de Phi no elemento padrao
    // apos o calculo do Gradiente
    // *************************************************
    // Aloca memoria para temp
    
    double grad[ndim][NGQP];
    double *ptr_grad[ndim];
    for(int k=0;k<ndim;++k){
      ptr_grad[k]=grad[k];
    }
    // Calcula o Laplaciano
    for(j=0;j<nn;++j){
      for(int l=0;l<NGQP;++l) LaplacianoPhi[i][j][l] = 0.0;
      for(int k=0;k<ndim;k++) { // direcao k
        // Calcular as derivadas de GradPhi[i][j][k]
        ptr_stdel[i]->Gradiente(ptr_grad,GradPhi[i][j][k],ptvert,Vert_map);
        for(int l=0;l<NGQP;l++) LaplacianoPhi[i][j][l] += grad[k][l];
      }
    }
    
    // *********************************************
    // Calcula os tracos de Phi e de seu gradiente
    // no elemento padrao
    // *********************************************
    ptr_stdel[i]->elem_traces(ptvert,Vert_map,sinal,TP,TrGradPhi[i],Jb);
    
    // Implementar nos elementos padroes (Triangle, Linear, LinearLeg, Quadrilateral e Tetrahedral)
    
    for(h=0;h<numborders;h++){
      for(int ndir=0;ndir<ndim;ndir++){
        n_e[ndir]=border[border_num[h]].normal[ndir];
      }
      for(j=0;j<nn;j++){
       	if(ptr_stdel[i]->is_on_border(j,h,pos)){
          for(int q=0;q<qmax;q++){
            TrPhi[i][h][pos][q] = TP[h][j][q];
          }
        }
        
        //	int ini = h*s0 + j*s1;
        for(int q=0;q<qmax;q++){
          // printf("q %3d/%3d ",q,qmax);
          double temp = 0.0;
          for(int ndir=0; ndir < ndim;ndir++){
            //	    TrGradPhi[i][h][j][ndir][q] = TGP[ini+ndir*qmax+q];
            temp+= (perm[ndir]*n_e[ndir])*TrGradPhi[i][h][j][ndir][q]; // gphi_[ndir][q]);
            // printf("Concluido TrGradPhi[%d][%d][%d][%d][%d] = % g\n",i,h,j,ndir,q,TrGradPhi[i][h][j][ndir][q]);
          }
          TrKgradPhi_n[i][h][j][q]=temp;
          //printf("\nTrKgradPhi_n[%d][%d][%d][%d] = % g\n\n",i,h,j,q,TrKgradPhi_n[i][h][j][q]);
        }
      }
    }
    
    // Libera memoria dinamica de TP
    for(h=0;h<numborders; h++) {
      for(j=0;j<nn;j++){
        delete [] TP[h][j]; TP[h][j]=nullptr;
      }
      delete [] TP[h]; TP[h] = nullptr;
    }
    delete [] TP; TP =nullptr;
  }
};
// ***************************************************************************
void DG_Elem::set_porosidade(double value)
{
  porosidade=value;
};
// ***************************************************************************
// ***************************************************************************
void DG_Elem::inicia_tracos( EDGE * border)
{
  // OBSOLETA!!!
  // Esta procedure necessita que haja pontos de Gauss nas bordas.
  // Para pontos de Gauss totalmente internos ou casos gerais usar
  // inicia_funcoes_na_borda (ver acima)
  // ******************************
  // Calcular Phi e seu gradiente *
  // ******************************
  int i,j,h,pos;
  double n_e[ndim];
  
  for(i=0;i<NumLocalVars;i++) {//i= variavel
    int nn=numn[i];
    int NGQP=ptr_stdel[i]->NGQP_val();
    int qmax = ptr_stdel[i]->qborder_val();
    double phi[NGQP];
    double gphi_[ndim][qmax];
    
    for(j=0; j < nn; j++) {// j=modo
      // calcular os valores de Phi_j nos pontos de quadratura
      ptr_stdel[i]->eval_Phi(j,phi);
      // Calcular as derivadas de Phi_j
      ptr_stdel[i]->Gradiente(GradPhi[i][j],phi,ptvert,Vert_map);
      
      // calcula os tracos de gradiente de Phi
      for(h=0; h < numborders; h++) { // h = aresta
        //printf("Ponto teste em inicia_tracos i=%d j=%d h=%d \n",i,j,h);
        for(int ndir=0; ndir<ndim;ndir++){
          // trace eh válido somente quando ha pontos de Gauss na borda
          ptr_stdel[i]->trace(h,qmax,sinal[h],GradPhi[i][j][ndir],gphi_[ndir]);
          n_e[ndir]=border[border_num[h]].normal[ndir];
        }
        for(int q=0; q < qmax; q++){
          // printf("q %3d/%3d ",q,qmax);
          double temp = 0.0;
          for(int ndir=0; ndir < ndim;ndir++){
            temp += (perm[ndir]*n_e[ndir]*gphi_[ndir][q]);
            TrGradPhi[i][h][j][ndir][q]=gphi_[ndir][q];
          }
          TrKgradPhi_n[i][h][j][q]=temp;
          //printf("Concluido\n");
        }
        // armazena o traco de phi
        if(ptr_stdel[i]->is_on_border(j,h,pos)){
          // printf("armazena traco de phi i=%d j= %d h=%d pos= %d ...\n",i,j,h,pos);
          // trace eh válido somente quando ha pontos de Gauss na borda
          ptr_stdel[i]->trace(h,qmax,sinal[h],phi,TrPhi[i][h][pos]);
          
          //printf(" salvou\n ");
          
        }
        //printf("Concluido h =%d (ne)\n",h);
      }
      //printf(" Concluido j = %d \n",j);
    }
    //printf(" Concluido i = %d \n",i);
  }
};

// ****************************************************************
// Calcula os tracos de funcoes somando os tracos dos modos
// O esforco computacional eh proporcional ao numero de modos
// ****************************************************************
void DG_Elem::calcula_tracos(Fluids fls)
{
  
  const int qborder=ptr_stdel[0]->qborder_val();
  
  const int ns=numn[sat];
  const int np=numn[pres];
  
  int pos;
  double aux;
  
  // Trsn
  for(int h=0; h < numborders; h++){
    for(int q=0; q < qborder; q++){
      aux=0.0;
      for(int n=0; n < ns; n++){
        if(ptr_stdel[sat]->is_on_border(n,h,pos)){
          aux+= TrPhi[sat][h][pos][q]*u0[sat][n];
        }
      }
      Trsn[h][q]=aux;
    }
  }
  
  // Trpw
  for(int h=0; h < numborders; h++){
    for(int q=0; q < qborder; q++){
      aux=0.0;
      for(int n=0; n < np ;n++){
        if(ptr_stdel[pres]->is_on_border(n,h,pos)){
          aux+= TrPhi[pres][h][pos][q]*u0[pres][n];
        }
      }
      Trpw[h][q]=aux;
    }
  }
  
  // TrKgrad_sn e TrKgrad_pc
  
  for(int h=0; h < numborders; h++){
    for(int q=0; q < qborder; q++){
      double aux1= fls.dpc(Trsn[h][q]);   // <---- fonte de dificuldade
      for(int ndir=0; ndir < ndim; ndir++){
        aux=0.0;
        for(int n=0; n < ns; n++){
          aux+=  TrGradPhi[sat][h][n][ndir][q]*u0[sat][n];
        }
        TrKgrad_sn[h][ndir][q]=aux*perm[ndir];
        TrKgrad_pc[h][ndir][q]=aux1*aux;
      }
    }
  }
  
  // TrKgrad_pw
  
  for(int h=0; h < numborders; h++){
    for(int q=0; q < qborder; q++){
      for(int ndir=0; ndir < ndim; ndir++){
        aux=0.0;
        for(int n=0; n < np; n++){
          aux += TrGradPhi[pres][h][n][ndir][q]*u0[pres][n];
        }
        TrKgrad_pw[h][ndir][q]=aux*perm[ndir];
      }
    }
  }
  
  /*
   //2
   
   // *************************************
   // Passar os tracos de sn e pw para os
   // vetores globais
   // *************************************
   if(gbtrpw!=NULL && gbtrsn!= NULL){
   int ind = gbtrbmap;
   for(int h=0;h<nborder;++h) {
			for(int q=0;q<qborder;++q) {
   gbtrpw[ind]=Trpw[h][q];
   gbtrsn[ind]=Trsn[h][q];
   ind++;
			}
   }
   }
   */
};
// *********************************************************************************
void DG_Elem::echo_traco(FILE * f_eco)
{
  int pos;
  for(int i=0;i<NumLocalVars;++i){
    int nn= numn[i];
    int qmax =  ptr_stdel[i]->qborder_val();
    
    for(int h=0;h<numborders;h++){
      
      for(int j=0;j<nn;j++){
        
        if(ptr_stdel[i]->is_on_border(j,h,pos)){
          for(int q=0;q<qmax;++q){
            printf("\nTrPhi[%d][%d][%d(%d)][%d] = % g\n\n",i,h,pos,j,q,TrPhi[i][h][pos][q]);
          }
        }
        for(int q=0;q<qmax;++q){
          // printf("q %3d/%3d ",q,qmax);
          // double temp = 0.0;
          for(int ndir=0; ndir < ndim;++ndir){
            printf("Concluido TrGradPhi[%d][%d][%d][%d][%d] = % g\n",i,h,j,ndir,q,TrGradPhi[i][h][j][ndir][q]);
          }
          printf("\nTrKgradPhi_n[%d][%d][%d][%d] = % g\n\n",i,h,j,q,TrKgradPhi_n[i][h][j][q]);
        }
      }
    }
  }
  
  if(f_eco != NULL){
    // *************************************************************
    // Echo dos resultados
    // *************************************************************
    int NGQP=ptr_stdel[0]->NGQP_val();
    fprintf(f_eco,"Jacobiano\n");
    for(int i=0;i<NGQP;++i)fprintf(f_eco,"%11.4e\n",JV[i]);
    
    fprintf(f_eco,"\nPhi para o elemento\n");
    for(int i=0;i<NumLocalVars;++i){
      fprintf(f_eco,"\n variavel %d\n",i);
      int nn     =numn[i];
      NGQP       =ptr_stdel[i]->NGQP_val();
      int qmax   =ptr_stdel[i]->qborder_val();
      
      for(int h=0;h<numborders;++h){
        fprintf(f_eco,"\naresta local %d sinal = %d\n",h,sinal[h]);
        for(int j=0;j<nn;++j){
          
          if(ptr_stdel[i]->is_on_border(j,h,pos)){
            fprintf(f_eco,"modo %d: ",j);
            for(int q=0;q<qmax;++q){
              //fprintf(f_eco,"\nTrPhi[%d][%d][%d(%d)][%d] = % g\n\n",i,h,pos,j,q,TrPhi[i][h][pos][q]);
              fprintf(f_eco,"%g ",TrPhi[i][h][pos][q]);
            }
            fprintf(f_eco,"\n");
          }
          /*
           for(int q=0;q<qmax;++q){
           // printf("q %3d/%3d ",q,qmax);
           // double temp = 0.0;
           for(int ndir=0; ndir < ndim;++ndir){
           fprintf(f_eco,"Concluido TrGradPhi[%d][%d][%d][%d][%d] = % g\n",i,h,j,ndir,q,TrGradPhi[i][h][j][ndir][q]);
           }
           fprintf(f_eco,"\nTrKgradPhi_n[%d][%d][%d][%d] = % g\n\n",i,h,j,q,TrKgradPhi_n[i][h][j][q]);
           }
           */
        }
      }
      /*
       for(int j=0;j<nn;++j){
       fprintf(f_eco,"         modo %d\n",j);
       fprintf(f_eco,"Grad_Phi[%d][%d]\n",i,j);
       for(int m=0;m<NGQP;++m){
       fprintf(f_eco," m= %3d %11.4e  %11.4e\n",m,GradPhi[i][j][0][m],GradPhi[i][j][1][m]);
       }
       //printf("AQUI1 ***************************************************\n");
       for(int h=0;h<nborder;++h){
       fprintf(f_eco,"  Traco na aresta %d \n",h);
       for(int m=0;m<qmax;++m){
       fprintf(f_eco,"%11.4e %11.4e\n",TrGradPhi[i][h][j][0][m],TrGradPhi[i][h][j][1][m]);
       }
       }
       //printf("AQUI2 **************************************************\n");
       }
       */
    }
    // **********************************************************
    // Fim do echo dos tracos de Phi e GradPhi
    // **********************************************************
  }
  
};
// ****************************************************************************
void DG_Elem::perm_val(double K[3])
{
  for(int i=0;i<3;i++)K[i]=perm[i];
};

// ****************************************************************************
// ****************************************************************************
void DG_Elem::Traco_sn(const int & h,double * saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  for(int q=0;q<qmax;q++)saida[q]=Trsn[h][q];
};

// ****************************************************************************
void DG_Elem::Traco_pw(const int & h,double *saida)
{
  int qmax=ptr_stdel[1]->qborder_val();
  for(int q=0;q<qmax;q++)saida[q]=Trpw[h][q];
};

// ****************************************************************************
void DG_Elem::Traco_Kgrad_pw(const int & h,double ** saida)
{
  int qmax=ptr_stdel[1]->qborder_val();
  for(int k=0;k<ndim;k++){
    for(int q=0;q < qmax;q++){
      saida[k][q]=TrKgrad_pw[h][k][q];
    }
  }
};
// ****************************************************************************
void DG_Elem::Traco_Kgrad_pc(const int & h,double ** saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  for(int k=0;k<ndim;k++){
    for(int q=0;q < qmax;q++) {
      saida[k][q]=TrKgrad_pc[h][k][q];
    }
  }
};
// ****************************************************************************
void DG_Elem::Traco_Kgrad_sn(const int & h,double ** saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  
  for(int k=0;k < ndim;k++){
    for(int q=0;q < qmax;q++){
      if(std::isnan(TrKgrad_sn[h][k][q]))cout << "Float was Not a Number: TrKgrad " << q << endl;
      saida[k][q]=TrKgrad_sn[h][k][q];
    }
  }
};
// ****************************************************************************
void DG_Elem::Traco_phi(const int & lado,const int & ivar,
                                     const int & ind, double * saida)
{
  int qmax=ptr_stdel[ivar]->qborder_val();
  int pos;
 
  if(ptr_stdel[ivar]->is_on_border(ind,lado,pos)){
    for(int q=0;q<qmax;q++){
      saida[q]=TrPhi[ivar][lado][pos][q];
    }
  }
  else {
    for(int q=0;q<qmax;q++){
      saida[q]=0.0;
    }
  }
};
// ****************************************************************************
void DG_Elem::Traco_grad_phi(const int & lado,const int & ivar,
                                          const int & ind,double ** saida)
{// Ficou desnecessaria quando criei Traco_Kgrad_phi_n abaixo
  int qmax=ptr_stdel[0]->qborder_val();
  for(int k=0;k<ndim;k++){
    for(int q=0;q<qmax;q++)
      saida[k][q]=TrGradPhi[ivar][lado][ind][k][q];
  }
};

// ****************************************************************************
void DG_Elem::Traco_Kgrad_phi_n(const int & lado,const int & ivar,
                                             const int & ind,double * saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  for(int q=0;q<qmax;q++)
    saida[q]=TrKgradPhi_n[ivar][lado][ind][q];
};
// ****************************************************************************
// ****************************************************************************
void DG_Elem::VolumeIntegrals_IG(Fluids fls,  int & count,
                                              int * Ti,int * Tj, double * Tx,
                                              double * B)
{
  //cout << "VolumeIntegrals_IG  \n";
  
  int qmax=ptr_stdel[sat]->qborder_val();
  int NGQP=ptr_stdel[sat]->NGQP_val();
  
  
  int np;//ns
  int h,k,m;
  double aux;
  double res0[NGQP];
  double res[ndim][NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  for(k=0;k<ndim;k++){
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP){
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  for(m=0;m<NGQP;m++)JW[m]*=JV[m];
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++)pc[m]=fls.pressao_capilar(sn[m]);
  ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
  
  for(m=0;m<NGQP;m++){
    // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
    for(k=0;k<ndim;k++){// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
    }
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
  }
  
  
  //printf("VolumeIntegrals_IG Ate aqui 1 ...\n");
  
  /*
   // salva os tracos de K * gradientes
   calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
   double test[nborder][ndim][qmax];
   for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   test[h][k][q]=TrKgrad_sn[h][k][q];
   }
   }
   }
   */
  
  for(h=0;h<numborders;h++){
    ptr_stdel[sat] ->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++){
      ptr_stdel[sat] ->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat] ->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  
  //
  // printf("VolumeIntegrals_IG Ate aqui 2\n");
  np=numn[pres];
  
  for(int rp=0;rp<np;rp++){ // PE_VI
    // *****************************
    // calculo do vetor
    // *****************************
    
    // integral 2 com sinal trocado
    for(int k=0; k < ndim; k++) {
      prodvv(NGQP,Kgrad_pc[k],GradPhi[pres][rp][k],res[k]);
    }
    for(int k=1; k < ndim; k++) {
      somavv(NGQP,res[0],res[k],res[0]);
    }
    prodvv(NGQP,lambdan,res[0],res0);
    integral(NGQP,JW,res0,aux);
    if(std::isnan(aux))cout << "VolumeIntegrals_IG  Float was Not a Number: aux " << aux << endl;
    B[gbnmap[pres][rp]] -= aux; // primeiro valor de B
    
    // Fim do vetor em PE_VI
    
    // *****************************
    // Matriz
    // *****************************
    for(int lp=0;lp<np;lp++){ // PE_VI_PP
      for(int k=0; k < ndim; k++) {
        prodev(NGQP,perm[k],GradPhi[pres][lp][k],res[k]);
      }
      for(int k=0; k < ndim; k++) {
        prodvv(NGQP,res[k],GradPhi[pres][rp][k],res[k]);
      }
      for(int k = 1; k < ndim; k++) {
        somavv(NGQP,res[0],res[k],res[0]);
      }
      prodvv(NGQP,lambdat,res[0],res0);
      integral(NGQP,JW,res0,aux);
      if(std::isnan(aux))cout << "VolumeIntegrals_IG  Float was Not a Number 2: aux " << aux << endl;
      Ti[count]=gbnmap[pres][rp];
      Tj[count]=gbnmap[pres][lp];
      Tx[count]=aux;
      count++;
    }
  }
  // OK 22/05/2008
};
// ****************************************************************************
// Volume Integrals for regular elements creating UMFPACK arrays
// ****************************************************************************
void DG_Elem::VolumeIntegrals_UMFPACK(const double Dt,Fluids fls, int & count,
                                                   int * Ti,int * Tj, double * Tx,
                                                   double * B,
                                                   double * gbtrsn, double * gbtrpw)
{
  //cout << "DG_VI.cc DG_Elem::VolumeIntegrals_UMFPACK\n";
  
  const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  
  int h,k,m;
  double aux,aaux;
  double res0[NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  double phi_r[NGQP],phi_l[NGQP]; // em PhElem
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  
  
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  
  for(k=0;k<ndim;k++) {
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  // calcula os pesos de Gauss no vetor JW
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  
  // multiplica os pesos de Gauss pelo Jacobiano JV
  for(m=0;m<NGQP;m++){
    JW[m]*=JV[m];
    // printf(" JW  = %g    JV = %g \n", JW[m], JV[m]);
  }
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++) {
    pc[m]=fls.pressao_capilar(sn[m]);
  }
  ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
  
  // ***********************************************************************
  // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
  for(m=0;m<NGQP;m++) {
    for(k=0;k<ndim;k++) {// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
    }
    
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
  }
  
  // *********************************************************************
  // salva os tracos de K * gradientes
  
  
  /*
   //calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
   double test[nborder][ndim][qmax];
   for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   test[h][k][q]=TrKgrad_pw[h][k][q];
   }
   }
   }
   */
  for(h=0;h<numborders;h++) {
    ptr_stdel[sat]->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++) {
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  
  // *************************************
  // Passar os tracos de sn e pw para os
  // vetores globais
  // *************************************
  int ind = gbtrbmap;
  for(h=0;h<numborders;h++) {
    for(int q=0;q<qmax;q++) {
      gbtrpw[ind]=Trpw[h][q];
      gbtrsn[ind]=Trsn[h][q];
      ind++;
    }
  }
  
  const int ns=numn[sat];
  const int np=numn[pres];
  const int ntot=ns+np;// Num total de modos no elemento
  const int n0 = numn[0];// Num de modos para a variavel 0
  
  // vetor com indices
  int mr,ml;
  double  mx [ntot*ntot];
  for(int i=0;i<ntot*ntot;i++) {
    mx[i]=0.0;
  }
  
  for(int rp=0;rp<np;rp++) { // PE_VI
    mr=indice(n0,pres,rp);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    aux=0.0;
    double btemp;
    for(int q=0;q<NGQP;q++) { // PE_VI
      
      double temp_pw = 0.0;
      double temp_pc = 0.0;
      for (int k = 0; k < ndim; k++) {
        temp_pw += Kgrad_pw[k][q]*GradPhi[pres][rp][k][q];
        temp_pc += Kgrad_pc[k][q]*GradPhi[pres][rp][k][q];
      }
      
      aux+=JW[q]*(lambdat[q]*temp_pw + lambdan[q]*temp_pc);
    }
    
    btemp=aux; // inicializa valor de B em PE_VI
    
    // termo de fonte
    aaux=qn+qw;
    if(aaux!=0.0) {
      ptr_stdel[pres]->eval_Phi(rp,phi_r);
      aux=0.0;
      for(int q=0;q<NGQP;q++)
        aux+=JW[q]*phi_r[q];
      //B[gbnmap[pres][rp]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    B[gbnmap[pres][rp]]+=btemp;
    // Fim do vetor em PE_VI
    
    // *****************************
    // Matriz
    // *****************************
    for(int lp=0;lp<np;lp++) { // PE_VI_PP
      ml=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double temp = 0.0;
        for (int k = 0; k < ndim; k++) {
          temp += perm[k]*GradPhi[pres][lp][k][m]*GradPhi[pres][rp][k][m];
        }
        
        aux+=( temp * lambdat[m]*JW[m]);
      }
      mx[npos(ntot,mr,ml)]+=aux;
    } // Fim de PE_VI_PP
    
    for(int ls=0;ls<ns;ls++) { // PE_VI_SP
      ml=indice(n0,sat,ls);
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double Vaux[ndim];
        for(k=0;k<ndim;k++) {
          Vaux[k]=( d_lambdat[m]*Kgrad_pw[k][m]+
                   d_lambdan[m]*Kgrad_pc[k][m]+
                   lambdan[m]*d2_pc[m]*Kgrad_sn[k][m] ) *
          phi_l[m]+
          lambdan[m]*d_pc[m]*perm[k]*GradPhi[sat][ls][k][m];
        }
        
        aaux=0.0;
        for(k=0;k<ndim;k++) {
          aaux+=Vaux[k]*GradPhi[pres][rp][k][m];
        }
        
        aux+=(aaux*JW[m]);
      }
      mx[npos(ntot,mr,ml)]+=aux;
    }// Fim de PE_VI_SP
  } // Fim de PE_VI
  
  for(int rs=0;rs<ns;rs++) { // SE_VI
    mr=indice(n0,sat,rs);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    double btemp;
    ptr_stdel[sat]->eval_Phi(rs,phi_r);
    
    aux=0.0;
    for(m=0;m<NGQP;m++) {
      aux+=JW[m]*(sn[m]-sna[m])*phi_r[m];
    }
    aux*=(-porosidade/Dt);
    
    // segunda integral
    aaux=0.0;
    for(m=0;m<NGQP;m++) {
      
      double temp = 0.0;
      for (int k =0; k <ndim; k++) {
        temp += Kgrad_pw[k][m]*GradPhi[sat][rs][k][m];
      }
      aaux+=JW[m]*lambdaw[m]*temp;
    }
    
    aux+=aaux;
    //B[gbnmap[sat][rs]]+=aux; // inicializa valor de B em SE_VI
    btemp=aux;
    // termo de fonte
    
    aaux=qw;
    
    if(aaux!=0.0) {
      integral(NGQP,JW,phi_r,aux);
      //B[gbnmap[sat][rs]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    B[gbnmap[sat][rs]] += btemp ;
    // *************************
    // Matriz
    // *************************
    
    for(int lp=0;lp<np;lp++) { // SE_VI_PS
      ml=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double temp = 0.0;
        for (int k = 0; k < ndim; k++) {
          temp += perm[k]*GradPhi[pres][lp][k][m]*GradPhi[sat][rs][k][m];
        }
        
        aux+=(temp * lambdaw[m]*JW[m]);
      }
      mx[npos(ntot,mr,ml)]+=aux;
    }
    
    for(int ls=0;ls<ns;ls++) { // SE_VI_SS
      ml=indice(n0,sat,ls);
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        aux+=JW[m]*(d_lambdaw[m]*res0[m]-porosidade/Dt*phi_r[m])*phi_l[m];
      }
      mx[npos(ntot,mr,ml)]+=aux;
    }
  }
  // ***************************************************************************
  // Passar a matriz mx para os vetores Ti, Tj e Tx
  // ***************************************************************************
  for(int rp=0;rp<np;rp++) {
    mr=indice(n0,pres,rp);
    
    for(int lp=0;lp<np;lp++) {
      ml=indice(n0,pres,lp);
      Ti[count]=gbnmap[pres][rp];
      Tj[count]=gbnmap[pres][lp];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
    
    for(int ls=0;ls<ns;ls++) {
      ml=indice(n0,sat,ls);
      Ti[count]=gbnmap[pres][rp];
      Tj[count]=gbnmap[sat][ls];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
  }
  for(int rs=0;rs<ns;rs++) {
    mr=indice(n0,sat,rs);
    
    for(int lp=0;lp<np;lp++) {
      ml=indice(n0,pres,lp);
      Ti[count]=gbnmap[sat][rs];
      Tj[count]=gbnmap[pres][lp];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
    
    for(int ls=0;ls<ns;ls++) {
      ml=indice(n0,sat,ls);
      Ti[count]=gbnmap[sat][rs];
      Tj[count]=gbnmap[sat][ls];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
  }
  // delete [] mx; mx=nullptr;
};
// alterado em 21/09/2011
// *****************************************************************************

// *****************************************************************************
// Usar se precisar atualizar so o termo de DT nas integrais de volume
//	 da equacao de sn
// *****************************************************************************
void DG_Elem::VolumeIntegralsT(const double Dt_new, const double Dt_old,
                                            int & count,
                                            int * Ti,int * Tj, double * Tx,
                                            double * B)
{
  //cout << "DG_VI.cc  DG_Elem::VolumeIntegralsT\n";
  
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  
  int m;
  double aux;
  double JW[NGQP];
  double sn[NGQP];
  double phi_r[NGQP]; // em PhElem
  
  const double fator=(1/Dt_new - 1/Dt_old);  // Para simples correcao do termo contendo Dt
  
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  
  for(m=0;m<NGQP;m++)JW[m]*=JV[m];
  
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  
  const int ns=numn[sat];
  
  for(int rs=0;rs<ns;rs++) { // SE_VI
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    
    ptr_stdel[sat]->eval_Phi(rs,phi_r);
    
    aux=0.0;
    for(m=0;m<NGQP;m++) {
      aux+=JW[m]*(sn[m]-sna[m])*phi_r[m];
    }
    aux*=(-porosidade*fator);
    
    B[gbnmap[sat][rs]]+=aux; // corrige valor de B em SE_VI
    
    
    // *************************
    // Matriz
    // *************************
    
    for(int ls=0;ls<ns;ls++) { // SE_VI_SS
      Ti[count]=gbnmap[sat][rs];
      Tj[count]=gbnmap[sat][ls];
      Tx[count]=(-porosidade*fator)*Mass_sn[rs][ls];//aux;
      count++;
    }
  }
};

// ****************************************************************************
// Volume Integrals for ghost elements - calcula os tracos
// ****************************************************************************
void DG_Elem::VolumeTracos(const double Dt,Fluids fls, double * gbtrsn, double * gbtrpw)
{
  //cout << "DG_VI.cc DG_Elem::VolumeTacos\n";
  
  const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  
  int h,k,m;
  double aux;
  //  double res0[NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  // double phi_r[NGQP],phi_l[NGQP];
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  
  for(k=0;k<ndim;k++) {
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  for(m=0;m<NGQP;m++) JW[m]*=JV[m];
  
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++)pc[m]=fls.pressao_capilar(sn[m]);
  ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
  // ***********************************************************************
  // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
  for(m=0;m<NGQP;m++) {
    for(k=0;k<ndim;k++) {// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
    }
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
  }
  // *********************************************************************
  // salva os tracos de K * gradientes
  //calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
  
  for(h=0;h<numborders;h++) {
    ptr_stdel[sat]->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++) {
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  
  //2
  
  // *************************************
  // Passar os tracos de sn e pw para os
  // vetores globais
  // *************************************
  int ind = gbtrbmap;
  for(h=0;h<numborders;h++) {
    for(int q=0;q<qmax;q++) {
      gbtrpw[ind]=Trpw[h][q];
      gbtrsn[ind]=Trsn[h][q];
      ind++;
    }
  }
};

/******************************************************************************/
/*  Uso de Trilinos      */
/******************************************************************************/
// ****************************************************************************
// Volume Integrals for regular elements creating Epetra_FECrsMatrix and
// Epetra_FEVector
// ****************************************************************************
void DG_Elem::VolumeIntegrals(const double Dt,Fluids fls,
                              Teuchos::RCP<Epetra_FECrsMatrix>  A,
                              Teuchos::RCP<Epetra_FEVector> RHS,
                              double * gbtrsn, double * gbtrpw)
{
  //cout<< "DG_VI.cc  DG_Elem::VolumeIntegrals\n";
  //#define indice(n0,var,i) ((n0*var)+i)
  //#define npos(ntot,i,j) ((ntot*i)+j)
  
  
  const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  
  int h,k,m;
  double aux,aaux;
  double res0[NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  double phi_r[NGQP],phi_l[NGQP];
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  
  for(k=0;k<ndim;k++) {
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  for(m=0;m<NGQP;m++)
    JW[m]*=JV[m];
  
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++){
    pc[m]=fls.pressao_capilar(sn[m]);
    }
    ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);// <-- derivada por colocacao; ptr_gpc esta no espaco aproximante
  for(m=0;m<NGQP;m++) {
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
    for(k=0;k<ndim;k++) {// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
      // ***********************************************************************
      // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
      // ***********************************************************************
      
      // ---------------------------------------------------------------------------------------------------------
      // Essa opcao gera problemas de convergencia
      //Kgrad_pc[k][m]=d_pc[m]*Kgrad_sn[k][m]; // <----- 06/01/2014 \Nabla{pc} = dpc * \Nabla{sn}; nao esta no espaco aproximante
        }
    }
  // *********************************************************************
  // salva os tracos de K * gradientes
  // Usar uma das duas alternativas abaixo
  
  // *******************************************************************
  // Alternativa 1: Calcula os tracos somando os tracos das
  // funcoes elementares previamente calculados
  // O esforco computacional eh proporcional ao numero de modos e ao
  // numero de pontos de Gauss
  //calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
  // *******************************************************************
  /*
   double test[nborder][ndim][qmax];
   for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   test[h][k][q]=TrKgrad_pc[h][k][q];
   }
   }
   }
   */
  // ou
  
  // *******************************************************************
  // Alternativa 2:
  // Interpola as funcoes nos pontos de Gauss nas bordas
  // O esforco computacional eh proporcional ao numero de pontos de Gauss
  // e nao depende do numero de modos
  // *******************************************************************
  
  for(h=0;h<numborders;h++) {
    ptr_stdel[sat] ->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++) {
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);  // <-- tem que ser calculado assim
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  
  // Teste
 
  // Fim do calculo dos tracos
  // ****************************************************************************
  // 3
  
  // *************************************
  // Passar os tracos de sn e pw para os
  // vetores globais
  // *************************************
  // cout << "gbtrbmap = "<< gbtrbmap << "\n";
  int ind = gbtrbmap;
  for(h=0;h<numborders;h++) {
    for(int q=0;q<qmax;q++) {
      gbtrpw[ind]=Trpw[h][q];
      gbtrsn[ind]=Trsn[h][q];
      ind++;
    }
  }
  
  const int ns=numn[sat];
  const int np=numn[pres];
  const int ntot=ns+np;// Num total de modos no elemento
  const int n0 = numn[0];// Num de modos para a variavel 0
  
  int mi,mj;
  double mx /*= new double*/ [ntot*ntot];
  double B  /*= new double*/ [ntot];
  int  indx /*= new int*/ [ntot];
  for(int i=0;i<ntot*ntot;i++) {
    mx[i]=0.0;
  }
  for(int i=0;i<ntot;i++) {
    B[i]=0.0;
  }
  
  for(int rp=0;rp<np;rp++) { // PE_VI
    mi=indice(n0,pres,rp);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    aux=0.0;
    double btemp;
    for(int q=0;q<NGQP;q++) { // PE_VI
      
      double temp_pw = 0.0;
      double temp_pc = 0.0;
      
      for (int k = 0; k < ndim; k++) {
        temp_pw += Kgrad_pw[k][q]*GradPhi[pres][rp][k][q];
        temp_pc += Kgrad_pc[k][q]*GradPhi[pres][rp][k][q];
      }
      aux+=JW[q]*(lambdat[q]*temp_pw + lambdan[q]*temp_pc);
    }
    
    btemp=aux; // inicializa valor de B em PE_VI
    
    // termo de fonte
    aaux=qn+qw;
    if(aaux!=0.0) {
      ptr_stdel[pres]->eval_Phi(rp,phi_r);
      aux=0.0;
      for(int q=0;q<NGQP;q++)
        aux+=JW[q]*phi_r[q];
      //B[gbnmap[pres][rp]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    //B[gbnmap[pres][rp]]+=btemp;
    B[mi] += btemp;
    // Fim do vetor em PE_VI
    
    // *****************************
    // Matriz
    // *****************************
    for(int lp=0;lp<np;lp++) { // PE_VI_PP
      mj=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        double temp = 0.0;
        
        for(int k =0; k <ndim; k++) {
          temp += (perm[k]*GradPhi[pres][lp][k][m]*GradPhi[pres][rp][k][m]);
        }
        aux+=( temp * lambdat[m]*JW[m] );
      }
      mx[npos(ntot,mi,mj)]+=aux;
    } // Fim de PE_VI_PP
    
    for(int ls=0;ls<ns;ls++) { // PE_VI_SP
      mj=indice(n0,sat,ls);
      
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      //cout<< "   <- Volume Integrals Ia " << npos(ntot,mi,mj) << "\n";
      for(m=0;m<NGQP;m++) {
        
        double Vaux[ndim];
        for(k=0;k < ndim;k++) {
          //cout<< "   <- sat ls k m "  << sat << "," << ls<< "," << k<< ","<< m <<"\n";
          Vaux[k]=( d_lambdat[m]*Kgrad_pw[k][m]+
                   d_lambdan[m]*Kgrad_pc[k][m]+
                   lambdan[m]*d2_pc[m]*Kgrad_sn[k][m]
                   ) *
          phi_l[m] +
          lambdan[m]*d_pc[m]*perm[k]*GradPhi[sat][ls][k][m];
          //cout<< "Volume Integrals Ia k =" << k << "\n" ;
        }
        
        aaux=0.0;
        
        for(k=0;k<ndim;k++) {
          aaux+=Vaux[k]*GradPhi[pres][rp][k][m];
        }
        
        aux+=(aaux*JW[m]);
        
      }
      
      mx[npos(ntot,mi,mj)]+=aux;
      
    }// Fim de PE_VI_SP
  } // Fim de PE_VI
  
  for(int rs=0;rs<ns;rs++) { // SE_VI
    //   mi=indice(n0,sat,rs);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    double btemp;
    ptr_stdel[sat]->eval_Phi(rs,phi_r);
    
    aux=0.0;
    for(m=0;m<NGQP;m++) {
      aux+=JW[m]*(sn[m]-sna[m])*phi_r[m];
    }
    aux*=(-porosidade/Dt);
    
    // segunda integral
    aaux=0.0;
    for(m=0;m<NGQP;m++) {
      double temp = 0.0;
      for (k=0; k< ndim; k++) {
        temp += Kgrad_pw[k][m]*GradPhi[sat][rs][k][m];
      }
      aaux+=JW[m]*lambdaw[m]*temp;
      res0[m]=temp;
    }
    aux+=aaux;
    //B[gbnmap[sat][rs]]+=aux; // inicializa valor de B em SE_VI
    btemp=aux;
    // termo de fonte
    
    aaux=qw;
    if(aaux!=0.0) {
      integral(NGQP,JW,phi_r,aux);
      //B[gbnmap[sat][rs]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    //B[gbnmap[sat][rs]] += btemp ;
    mi = indice(n0,sat,rs);
    B[mi] += btemp ;
    // *************************
    // Matriz
    // *************************
    
    for(int lp=0;lp<np;lp++) { // SE_VI_PS
      mj=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double temp =0.0;
        for ( int k=0; k < ndim; k++) {
          temp += (perm[k]*GradPhi[pres][lp][k][m]*GradPhi[sat][rs][k][m]);
        }
        
        aux+=(temp *lambdaw[m]*JW[m]);
      }
      mx[npos(ntot,mi,mj)]+=aux;
    }
    
    for(int ls=0;ls<ns;ls++) { // SE_VI_SS
      mj=indice(n0,sat,ls);
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        aux+=JW[m]*(d_lambdaw[m]*res0[m]-porosidade/Dt*phi_r[m])*phi_l[m];
      }
      mx[npos(ntot,mi,mj)]+=aux;
    }
  }
  // ***************************************************************************
  // Colocar os indices no vetor Ti
  // ***************************************************************************
  for(int rp=0;rp<np;rp++) {
    mi=indice(n0,pres,rp);
    indx[mi] = gbnmap[pres][rp];
  }
  for(int rs=0;rs<ns;rs++) {
    mi=indice(n0,sat,rs);
    indx[mi] = gbnmap[sat][rs];
  }
  A->InsertGlobalValues(ntot,indx,mx,Epetra_FECrsMatrix::ROW_MAJOR);
  RHS->SumIntoGlobalValues(ntot,indx,B);
  //  for(int i = 0; i < ntot*ntot; i++) {
  //    printf("mx[%3d] = %12.5e\n",i,mx[i]);
  //  }
  
};
// alterado em 21/09/2011
// alterado em 21/10/2011
// alterado em 23/10/2011
// alterado em 13/02/2013

// **************************************************************
// Inicaliza os vetores locais: JV,b0,bs,Grad e traco
// **************************************************************

void DG_Elem::inicia_vetores()
{
  int nn;//nb,q0,q1;
  int h,i,j,k;
  
  if(vetores_iniciados == 1) { cout<< "vetores locais de PhElem já iniciados\n"; }
  
  else {
    vetores_iniciados = 1;
    //cout<< "DG_Elem::inicia_vetores()\n";
    for (k=0;k<NumLocalVars;++k){
      nn=numn[k];
      u0[k] = new double [nn];
      usave[k] = new double [nn];
      for(i=0;i<nn;++i){
        u0[k][i]=0.0;
      }
    }
    // **************************************
    // Calcular o Jacobiano   *             *
    // **************************************
    // q0= ptr_stdel[0]->Q_val(0);
    // q1= ptr_stdel[0]->Q_val(1);
    
    int qmax = ptr_stdel[0]->qborder_val();
    int NGQP;
    //printf("Calculo do Jacobiano: dimensao=%d\n",NGQP);
    //JV = new double [NGQP]; // opcao 2: um unico vetor
    //ptr_stdel[0]->Jacobian(ptvert,Vert_map,JV);
    compute_JV(0);
    
    sna = new double [ ptr_stdel[0]->NGQP_val() ];
    pwa = new double [ ptr_stdel[1]->NGQP_val() ];
    
    // ***********************************************************
    // Alocacao dinamica de memoria para a matriz de massa de sn *
    // Mass_sn[i][j]                                             *
    // ***********************************************************
    int ns=numn[0];
    Mass_sn = new double * [ns];
    for(i=0;i<ns;i++)
      Mass_sn[i] = new double [ns];
    // ******************
    // Calcular Mass_sn *
    // ******************
    for(i=0;i<ns;i++){
      Mass_sn[i][i]= ptr_stdel[0]->mass(i,i,JV);
      for(j=i+1;j<ns;j++){
        Mass_sn[i][j] = ptr_stdel[0]->mass(i,j,JV);
        Mass_sn[j][i]=Mass_sn[i][j];
      }
    }
    
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // do vetor Jacobiano Jb
    // [lado] [var] [indice] [direcao] [posicao]   *
    // *********************************************
    
    Jb = new double [numborders*qmax];
    //cout << "INICIOU JB !!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // **************************************
    // Alocacao dinamica de memoria para os *
    // gradientes dos modos                 *
    // [var] [modo] [direcao] [posicao]   *
    // **************************************
    GradPhi = new double *** [NumLocalVars];// variavel i NumLocalVars=2
    for(i=0;i<NumLocalVars;++i){
      nn=numn[i];
      NGQP=ptr_stdel[i]->NGQP_val();
      GradPhi[i] = new double ** [nn];// modo j
      for(j=0;j<nn;j++){
        GradPhi[i][j] = new double * [ndim];// direcao k
        for(k=0;k<ndim;k++){
          GradPhi[i][j][k] = new double [NGQP];// posicao
        }
      }
    }
    // **************************************
    // Alocacao dinamica de memoria para os *
    // Laplacianos  dos modos               *
    // [var] [modo] [posicao]             *
    // **************************************
    LaplacianoPhi = new double ** [NumLocalVars];// variavel i  NumLocalVars=2
    for(i=0;i<NumLocalVars;++i){
      nn  =numn[i];
      NGQP=ptr_stdel[i]->NGQP_val();
      LaplacianoPhi[i] = new double * [nn];// indice j
      for(j=0;j<nn;j++){
        LaplacianoPhi[i][j] = new double [NGQP];// posicao
      }
    }
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [var] [lado] [modo] [direcao] [posicao]   *
    // *********************************************
    TrGradPhi = new double **** [NumLocalVars];// variavel i  NumLocalVars=2
    for(i=0;i<NumLocalVars;++i){
      nn=numn[i];
      TrGradPhi[i] = new double *** [numborders];//lado h
      for(h=0;h<numborders;h++){
        TrGradPhi[i][h] = new double ** [nn];// modo j
        for(j=0;j<nn;j++){
          TrGradPhi[i][h][j] = new double * [ndim];// direcao k
          for(k=0;k<ndim;k++){
            TrGradPhi[i][h][j][k] = new double [qmax];// posicao
          }
        }
      }
    }
    
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [var] [lado] [modo] [posicao]               *
    // *********************************************
    TrKgradPhi_n = new double *** [NumLocalVars];//variavel i  NumLocalVars=2
    for(i=0;i<NumLocalVars;++i){
      nn=numn[i];
      TrKgradPhi_n[i] = new double ** [numborders];// lado h
      for(h=0;h<numborders;h++){
        TrKgradPhi_n[i][h] = new double * [nn];// modo j
        for(j=0;j<nn;j++){
          TrKgradPhi_n[i][h][j] = new double [qmax];// posicao
        }
      }
    }
    
    // ****************************************************
    // Alocacao dinamica de memoria para os tracos de phi *
    // [var] [lado] [posicao modo]  [posicao q]           *
    // ****************************************************
    TrPhi = new double *** [NumLocalVars];// variavel i  NumLocalVars=2
    for(i=0;i<NumLocalVars;++i){
      TrPhi[i] = new double ** [numborders]; //lado h
      for(h=0;h<numborders;h++){
        int pmax = 1;
        for(int t=0; t < (ndim-1); t++) {
          pmax *= ( (ptr_stdel[i]->P_val(t)) + 1);// pmax = P +1
        }
        TrPhi[i][h] = new double * [pmax];// posicao do modo j
        for(j=0; j < pmax; j++){
          TrPhi[i][h][j] = new double  [qmax];// posicao de q
        }
      }
    }
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [lado][posicao]                             *
    // *********************************************
    Trsn = new double * [numborders];//lado h
    for(h=0;h<numborders;h++){
      Trsn[h] = new double [qmax];// posicao
    }
    
    Trpw = new double * [numborders];//lado h
    for(h=0;h<numborders;h++){
      Trpw[h] = new double [qmax];// posicao
    }
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [lado] [direcao] [posicao]                  *
    // *********************************************
    TrKgrad_sn = new double ** [numborders];//lado h
    for(h=0;h<numborders;h++){
      TrKgrad_sn[h] = new double * [ndim];// direcao k
      for(k=0;k<ndim;k++){
        TrKgrad_sn[h][k] = new double  [qmax];// posicao
      }
    }
    TrKgrad_pw = new double ** [numborders];//lado h
    for(h=0;h<numborders;h++){
      TrKgrad_pw[h] = new double * [ndim];// direcao k
      for(k=0;k<ndim;k++){
        TrKgrad_pw[h][k] = new double  [qmax];// posicao
      }
    }
    TrKgrad_pc = new double ** [numborders];//lado h
    for(h=0;h<numborders;h++){
      TrKgrad_pc[h] = new double * [ndim];// direcao
      for(k=0;k<ndim;k++){
        TrKgrad_pc[h][k] = new double  [qmax];// posicao
      }
    }
  }
};

// ****************************************************************************
void DG_Elem::finaliza_vetores()
{
  if(vetores_iniciados == 1) {
    vetores_iniciados = 0;
    //cout<< "DG_Elem::finaliza_vetores()\n";
    delete [] JV; JV=nullptr;
    delete [] Jb; Jb=nullptr;
    
    int i,j,k;
    for (k=0;k<2;k++){
      delete [] u0[k]; u0[k]=nullptr;
      //   delete [] ua[k];ua[k]=nullptr;
      delete [] usave[k];usave[k]=nullptr;
    }
    delete [] stgbtrbmapM;stgbtrbmapM=nullptr;
    delete [] stgbtrbmapP;stgbtrbmapP=nullptr;
    
    delete [] sna;sna=nullptr;
    delete [] pwa;pwa=nullptr;
    
    int ns=numn[0];
    for(i=0;i<ns;i++){
      delete [] Mass_sn[i]; Mass_sn[i]=nullptr;
    }
    delete [] Mass_sn; Mass_sn=nullptr;
    
    for(i=0;i<2;i++){
      int nn=numn[i];
      for(j=0;j<nn;j++){
        for(k=0;k<ndim;k++){
          delete [] GradPhi[i][j][k]; GradPhi[i][j][k]=nullptr;// posicao
        }
        delete [] GradPhi[i][j];GradPhi[i][j]=nullptr;
        delete [] LaplacianoPhi[i][j]; LaplacianoPhi[i][j]=nullptr;
      }
      delete [] GradPhi[i];GradPhi[i]=nullptr;
      delete [] LaplacianoPhi[i]; LaplacianoPhi[i]=nullptr;
    }
    delete [] GradPhi; GradPhi=nullptr;
    delete [] LaplacianoPhi; LaplacianoPhi=nullptr;
    
    // traco de grad phi
    for(int i=0;i<2;i++){
      for(int h=0;h<numborders;h++){
        int nn=numn[i];
        for(j=0;j<nn;j++){
          for(k=0;k<ndim;k++){
            delete [] TrGradPhi[i][h][j][k];TrGradPhi[i][h][j][k]=nullptr;// posicao
          }
          delete [] TrGradPhi[i][h][j];TrGradPhi[i][h][j]=nullptr;
        }
        delete [] TrGradPhi[i][h];TrGradPhi[i][h]=nullptr;
      }
      delete [] TrGradPhi[i]; TrGradPhi[i]=nullptr;
    }
    delete [] TrGradPhi;TrGradPhi=nullptr;
    // traco de grad phi dot n
    for(int i=0;i<2;i++){
      for(int h=0;h<numborders;h++){
        int nn=numn[i];
        for(j=0;j<nn;j++){
          delete [] TrKgradPhi_n[i][h][j];TrKgradPhi_n[i][h][j]=nullptr;
        }
        delete [] TrKgradPhi_n[i][h];TrKgradPhi_n[i][h]=nullptr;
      }
      delete [] TrKgradPhi_n[i];TrKgradPhi_n[i]=nullptr;
    }
    delete [] TrKgradPhi_n;TrKgradPhi_n=nullptr;
    
    // traco de Phi
    for(int i=0;i < 2;i++){
      int p0=ptr_stdel[i]->P_val(0);
      for(int h=0;h < numborders;h++){
        for(j=0;j<p0;j++){
          delete [] TrPhi[i][h][j]; TrPhi[i][h][j]=nullptr;
        }
        delete [] TrPhi[i][h];TrPhi[i][h]=nullptr;
      }
      delete [] TrPhi[i];TrPhi[i]=nullptr;
    }
    delete [] TrPhi;TrPhi=nullptr;
    // *********************************************
    // Liberacao de memoria para os tracos         *
    // [lado][posicao]                             *
    // *********************************************
    int h;
    for(h=0;h<numborders;h++){
      delete [] Trsn[h];Trsn[h]=nullptr;
      delete [] Trpw[h];Trpw[h]=nullptr;
    }
    delete [] Trsn;Trsn=nullptr;
    delete [] Trpw;Trpw=nullptr;
    for(h=0;h<numborders;h++){
      for(k=0;k<ndim;k++){
        delete [] TrKgrad_sn[h][k];TrKgrad_sn[h][k]=nullptr;
        delete [] TrKgrad_pw[h][k];TrKgrad_pw[h][k]=nullptr;
        delete [] TrKgrad_pc[h][k];TrKgrad_pc[h][k]=nullptr;
      }
      delete [] TrKgrad_sn[h]; TrKgrad_sn[h]=nullptr;
      delete [] TrKgrad_pw[h]; TrKgrad_pw[h]=nullptr;
      delete [] TrKgrad_pc[h]; TrKgrad_pc[h]=nullptr;
    }
    delete [] TrKgrad_sn;TrKgrad_sn=nullptr;
    delete [] TrKgrad_pw;TrKgrad_pw=nullptr;
    delete [] TrKgrad_pc;TrKgrad_pc=nullptr;
  }
};
// ***************************************************************************************

