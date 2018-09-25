#include "spectral.h"
#include "DG_Prob.h"
// ****************************************************************************
void DG_Prob::DG_EI_flux(const EDGE border, double & flux_w, double & flux_n)
{
  //int m;
  MyElem e[2];
  const int iE=0;
  e[0]=el[border.elemento[0]];
  const int qmax=e[0].show_ptr_stdel(0)->qborder_val();// Pontos de gauss nas arestas
  const int ndim=e[0].show_ptr_stdel(0)->ndim_val();

  // substituir por uma funcao que retorne so os w's

  double length;
  double Jac[qmax],w[qmax],Dtemp[MAXQ][MAXQ];
  if(ndim==1){ 
    w[0]=1.0;
    length=1.0;
  }
  else {
    Gauss_Jacobi_parameters(qmax,0.0,0.0,Jac,w,Dtemp);
    length=border.modulo; 
  }
  // *******************************************************
  for(int q=0;q<qmax;q++) w[q] *= length/2.0; // w=w*length/2.0
  // *******************************************************
  double n_e[2];
  int a[2];
  double * sn[2];
  double lambdaw[2][qmax],lambdan[2][qmax],lambdat[2][qmax];
  double K_g_pw_n[ndim][qmax],K_g_pc_n[ndim][qmax];
  double mun=fluids.show_mun();
  double muw=fluids.show_muw();
  
  // double kn,kln,klw,klt;
  for(int i =0;i<ndim;i++) { 
    n_e[i]=border.normal[i];
  }
  
  a[iE]=border.num_local[0];
  
  // ************************************
  // Usa os tracos globais de sn e pw   *
  // ************************************
  int offset=border.gbtrbind[0];
  sn[0]=gbtrsn+offset;
  // pw[0]=gbtrpw+offset;
 
  for(int q = 0; q < qmax; q++){
    // Calculo dos escalares
    double aux=sn[iE][q];
    lambdan[iE][q]=fluids.Krn(aux)/mun;
    lambdaw[iE][q]=fluids.Krw(aux)/muw;
    lambdat[iE][q]=lambdan[iE][q]+lambdaw[iE][q];
  }

  e[iE].Traco_Kgrad_pw(a[iE],Kgpw[iE]);
  e[iE].Traco_Kgrad_pc(a[iE],Kgpc[iE]);
    
  for(int q=0;q<qmax;q++){
      K_g_pw_n[iE][q]=0.0;
      K_g_pc_n[iE][q]=0.0;
      for ( int k = 0; k < ndim; k++){ 
        K_g_pw_n[iE][q]+=(n_e[k]*Kgpw[iE][k][q]);
        K_g_pc_n[iE][q]+=(n_e[k]*Kgpc[iE][k][q]);
      }
  } 
  // ------------------------
  // Fluxos sobre as bordas
  // ------------------------
  flux_w=0.0;
  flux_n=0.0;
  for(int q = 0; q < qmax; q++){
    flux_w += w[q]*lambdaw[0][q]* K_g_pw_n[0][q];
    flux_n += w[q]*lambdan[0][q]*(K_g_pw_n[0][q]+K_g_pc_n[0][q]);
  }
};
