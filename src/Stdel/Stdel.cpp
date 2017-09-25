// $Id$

/**
 * @file
 * Stdel and Vertices class functions
 */

#include "spectral.h"
/*
// Class Vertice

Vertice::Vertice(double a,double b,double c)
{
  set_Vertice(a,b,c);
};

void Vertice::set_Vertice(double a,double b,double c)
{
  x=a;
  y=b;
  z=c;
};

double Vertice::x const
{
  return x;
};

double Vertice::y const
{
  return y;
};

double Vertice::z const
{
  return z;
};

void Vertice::print(FILE * file) const
{
  fprintf(file,"%11.4e %11.4e %11.4e\n",x,y,z);
}
*/

//  // Class Edge
//  Edge::Edge()
//  {
//    set(0,0,0);
//  };
//  
//  void Edge::set(int i, int j, int k)
//  {
//    Na=i;
//    Nb=j;
//    flag=k;
//  };

double  Stdel::show_Mb_comp(int i, int j) const
{
  return Mb_comp[i][j];
};

double  Stdel::show_Mi_inv(int i,int j) const
{
  return Mi_inv[i][j];
};

double  Stdel::show_McMi_inv(int i,int j) const
{
  return McMi_inv[i][j];
};

void Stdel::show_ind(int i,int& p,int& q,int& r)
{
  p=mode_[i].p_val();
  q=mode_[i].q_val();
  r=mode_[i].r_val();
}
//  
void Stdel::printStdel() const
{
  printf("ndim= %d\nP   Q  gqt\n",ndim);
  for (int i=0; i<ndim; i++)
    printf("%d   %d   %d\n",P[i], Q[i], gqt[i]);
};

void Stdel::print_matrices(FILE * fout)
{
  int i,j;
  fprintf(fout,"Stdel::print_matrices Matriz Mb_comp\n");  
  for(i=0;i<nb;i++){
    //fprintf(fout,"linha %d\n",i);
    for(j=0;j<nb;j++){
      fprintf(fout,"%11.4e ",Mb_comp[i][j]);
    }
    fprintf(fout,"\n");
  }
  
  fprintf(fout,"Stdel::print_matrices Matriz Mi_inv\n");  
  for(i=0;i<nn-nb;i++){
    //fprintf(fout,"linha %d\n",i);
    for(j=0;j<nn-nb;j++){
      fprintf(fout,"%11.4e ",Mi_inv[i][j]);
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"Stdel::print_matrices Matriz McMi_inv\n");  
  for(i=0;i<nb;i++){
    //fprintf(fout,"linha %d\n",i);
    for(j=0;j<nn-nb;j++){
      fprintf(fout,"%11.4e ",McMi_inv[i][j]);
    }
    fprintf(fout,"\n");
  }
};
// *********************************************************************
void Stdel::duplicar_mass_matrix(int NFields)
{
  int Nb=nb*NFields;
  int ni=(nn-nb);
  int Ni=ni*NFields;
  // Matrizes temporarias
  double Ab[Nb][Nb];
  double Ac[Nb][Ni];
  double Ai[Ni][Ni];
  // zerar as matrizes temporarias
  for(int i=0;i<Nb;i++){
    for(int j=0;j<Nb;j++)
      Ab[i][j]=0.0;
    for(int j=0;j<Ni;j++)
      Ac[i][j]=0.0;
  }
  for(int i=0;i<Ni;i++)
    for(int j=0;j<Ni;j++)
      Ai[i][j]=0.0;
  
  for(int a=0;a<NFields;a++){
    for(int i=0;i<nb;i++){
      int ii = i*NFields+a;
      for  (int j=0;j<nb;j++){
	int jj = j*NFields+a;
	Ab[ii][jj]=show_Mb_comp(i,j);
      }
      for(int j=0;j<ni;j++){
	int jj = j*NFields+a;
	Ac[ii][jj]=show_McMi_inv(i,j);
      }
    }
    for(int i=0;i<ni;i++){
      int ii = i*NFields+a;
      for(int j=0;j<ni;j++){
	int jj = j*NFields+a;
	Ai[ii][jj]=show_Mi_inv(i,j);
      }
    }
  }
  // Armazenar as matrizes Mb_comp, McMi_inv e Mi_inv
  for(int i=0;i<Nb;i++){
    for(int j=0;j<Nb;j++)
      Mb_comp[i][j]=Ab[i][j];
    for(int j=0;j<Ni;j++)
      McMi_inv[i][j]=Ac[i][j];
  }
  for(int i=0;i<Ni;i++)
    for(int j=0;j<Ni;j++)
      Mi_inv[i][j]=Ai[i][j];
};
// ****************************************************************************
int Stdel::w_vec(double saida[])
{
  int count=0;
  double aux2,aux1;
  for(int k=0;k<Q[2];k++){
    aux2=wGQ[2][k];
    for(int j=0;j<Q[1];j++){
      aux1=wGQ[1][j]*aux2;
      for(int i=0;i<Q[0];i++){
	saida[count++]=wGQ[0][i]*aux1; 
      }
    }
  }
  return count;
};



  
//  int Stdel::show_emapv(int i)
//  { 
//    return emapv[i];
//  };
//  int Stdel::show_emapi(int a)
//  { 
//    return emapi[a];
//  }
