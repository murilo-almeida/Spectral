
/*! \file ADifferentiation.cpp
  \brief  ADifferentiation.cpp documentação

  Este arquivo contém as rotinas para diferenciação 
  numérica, pelo método de Colocação de Gauss-Jacobi.

O Material aqui implementado está documentado no Apêndice C do livro
 "Spectral/hp Element Methods for CFD" escrito por
 G. Karniadakis & S. J. Spencer. 
*/

# include "spectral.h"
// ****************************************************************************
//! Collocation Differentiation
//! Gauss-Lobatto-Jacobi.
//! Implements material contained in
//! Appendix C of Spectral/hp Element Methods for CFD by G. Karniadakis & S. J. Spencer.
// ****************************************************************************
void ColDifGLJ(int Q, double alpha, double beta, double D[MAXQ][MAXQ])
{
  double dp1[MAXQ], dp2[MAXQ];
 // double sgn=1.0;
  double y, dy, d2y, aux, aaux;
  double x[Q];
  // Calcula os pontos de quadratura
  Jacobi_roots(Q-2,alpha+1.0, beta+1.0, x);
  for (int i=Q-2; i>0; i--) x[i]=x[i-1];
  x[0]   = -1.0;
  x[Q-1] =  1.0;

// ******************************************************************
// Forma alternativa
//  
//  gq1=gamma((double)Q-1.0);
//  gqa=gamma((double)Q+alpha);
//  gqb=gamma((double)Q+beta);
//  ga2=gamma(alpha+2.0);
//  gb2=gamma(beta+2.0);
//  if( (Q/2)*2 != Q)sgn=-1.0;
//  
//  dp1[0]=sgn*2.0*gqb/gq1/gb2;
//  dp2[0]=sgn*2*(alpha-((double)Q-1.0)*((double)Q+alpha+beta))/(beta+2.0)*gqb/gq1/gb2;
//  
//  dp1[Q-1]=-2.0*gqa/gq1/ga2;
//  dp2[Q-1]=-(beta-((double)Q-1)*((double)Q+alpha+beta))/(alpha+2.0)*dp1[Q-1];
//  for(int i=1; i<=Q-2; i++){
//    Jacobi_P(Q-2, alpha+1.0, beta+1.0, x[i], y, dy);
//    dp1[i]=(1.0-x[i])*(1+x[i])*dy;
//    dp2[i]=(alpha-beta+ (alpha+beta)*x[i]) * dy;
//  }
  
  for(int i=0; i<Q;i++){
    aux=x[i];
    aaux=1.0-(aux*aux);
    Jacobi_P(Q-2, alpha+1.0, beta+1.0, x[i], y, dy,d2y);
    dp1[i]=aaux*dy-2.0*aux*y;
    dp2[i]=aaux*d2y-4.0*aux*dy-2.0*y;
  }
  
  for(int i=0; i<Q; i++)
    for(int j=0; j<Q; j++){
      if(i==j) D[i][j] = dp2[i]/2.0/dp1[i];
      else
	D[i][j] = dp1[i]/dp1[j]/(x[i]-x[j]);
    }
}

// *******************************************************************
//
//! Collocation Differentiation
//! Gauss-Radau-Jacobi.
//! Implements material contained in
//! Appendix C of Spectral/hp Element Methods for CFD by G. Karniadakis & S. J. Spencer.
// *******************************************************************
void ColDifGRJ(int Q, double alpha, double beta, double D[MAXQ][MAXQ])
{
  double dp1[MAXQ], dp2[MAXQ];
  double y, dy, d2y, aux, aaux;
  double x[Q];
  // Calcula os pontos de quadratura
  Jacobi_roots(Q-1,alpha, beta+1.0, x);
  for (int i=Q-1; i>0; i--) x[i]=x[i-1];
  x[0]   = -1.0;
  
  for(int i=0; i<Q;i++){
    aux=x[i];
    aaux=1.0+aux;
    Jacobi_P(Q-1, alpha, beta+1.0, aux, y, dy,d2y);
    dp1[i]=aaux*dy+y;
    dp2[i]=aaux*d2y+2.0*dy;
  }
  
  for(int i=0; i<Q; i++)
    for(int j=0; j<Q; j++){
      if(i==j) D[i][j] = dp2[i]/2.0/dp1[i];
      else
	D[i][j] = dp1[i]/dp1[j]/(x[i]-x[j]);
    }
}

// *******************************************************************
//
//! Collocation Differentiation
//! Gauss-Jacobi.
//! Implements material contained in
//! Appendix C of Spectral/hp Element Methods for CFD by G. Karniadakis & S. J. Spencer.
// *******************************************************************

void ColDifGJ(int Q, double alpha, double beta, double D[MAXQ][MAXQ])
{
  double dp1[MAXQ], dp2[MAXQ];
//  double sgn=1.0;
  double y;
  double x[Q];
  // Calcula os pontos de quadratura
  Jacobi_roots(Q,alpha, beta, x);

  for(int i=0;i<Q;i++) Jacobi_P(Q,alpha,beta,x[i],y,dp1[i],dp2[i]);
  for(int i=0; i<Q; i++)
    for(int j=0; j<Q; j++){
      if(i==j) D[i][j] = dp2[i]/2.0/dp1[i];
      else
	D[i][j] = dp1[i]/dp1[j]/(x[i]-x[j]);
    }
}
