//
//  DG_EI_source.cpp
//  SDG
//
//  Created by Murilo Almeida on 11/30/13.
//  Copyright (c) 2013 Murilo Almeida. All rights reserved.
//
// DG (Discontinuois Galerkin) integrais de contorno (nos Edges)
// Calcula matriz mx e vetor B no contorno de entrada (Inflow)
//
#include "spectral.h"
#include "DG_Prob.h"

int DG_Prob::DG_EI_Inflow(const EDGE border,
													double * mx,
													double * B,
													int    * indx)
{
	
	// **************************************************************************
#include "DG_EI_Header.h" // Inicia variaveis locais comuns a DG_EI
	// **************************************************************************
	
  // **************************************************************************
  // Boundary Edges
  // **************************************************************************
	
  double * pdir, * sdir;
  pdir=border.pdir;
  //printf("Inside EdgeIntegrals: tipo = -1\n");
  sdir=border.sdir;
  double difpw[qmax],difpc[qmax];					        // PE V
  for (int q=0;q<qmax;q++){
    difpw[q]=pw[0][q]-pdir[q];
    difpc[q]=pc[0][q]-fluids.pressao_capilar(sdir[q]);
  }
	
	// **************************************************************************
	int I0,I1;
	iEr = 0;
	iEl = 0;
	const int ns=e[0].show_ptr_stdel(sat)->nn_val();// Num de modos para saturacao
  const int np=e[0].show_ptr_stdel(pres)->nn_val();// Num de modos para pressao
  const int ntot=ns+np;// Num total de modos no elemento
	int mi,mj;
  int iaux=ntot*ntot;
	
  for(int i=0; i<iaux; ++i) {
    mx[i]=0.0;
  }
  for(int i=0; i<ntot; ++i) {
    B[i]=0.0;
  }
	
	
  // *******************
  //     VETOR         *
  // *******************
  // *************************
  // Pressure Equation Vector
  // *************************
	
  //integral 1 do caderno
  pos=0;
  I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);		       // PE V
  I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);		       // PE V
  for(int ip=I0;ip<I1;ip++){					       // PE V
    int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip);
    e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);pos++;
    
    aux=0.0;
    for(int i=0;i<qmax;i++){
      aux+=w[i]*phi_r[i]*lambdat[0][i]*K_g_pw_n[0][i];
      if(std::isnan(phi_r[i]))cout << "Float was Not a Number: phi_r " << i<< endl;
      if(std::isnan(lambdat[0][i]))cout << "Float was Not a Number: lambdat " << i <<endl;
      if(std::isnan(K_g_pw_n[0][i]))cout << "Float was Not a Number: K_g_pw_n " << i<< endl;
    }
    mi=indice(n0,pres,rp);
    if(std::isnan(aux) )cout << "Float was Not a Number: ponto 1" << aux << endl;
		
    B[mi]  -= aux;
  }
  
  // integrais 3 e 7 do caderno (sao agrupadas numa so)
  
  for(int rp=0;rp<np;rp++){
    
    aux=0.0;
    for(int i=0;i<qmax;i++){
      aux+=w[i]*K_g_phi_n[iEr][pres][rp][i]*lambdat[iEr][i]*difpw[i];
    }
    
    mi=indice(n0,pres,rp);
    
    if(std::isnan(aux) )cout << "Float was Not a Number: ponto 2" << aux << endl;
    B[mi] += aux;
  }
  
  // integrais 5 e 8 do caderno  (sao agrupadas numa so)
  aux2=klt*penalty;
  if(std::isnan(aux2) )cout << "Float was Not a Number: ponto klt*penalty " << aux2 << " " << klt << " " << penalty << endl;
	
  I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);		       // PE V
  I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);
  //printf("I0 = %d    I1=%d\n",I0,I1);		       // PE V
  pos=0;
  for(int ip=I0;ip<I1;ip++){					       // PE V
    int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip);		       // PE V
    e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);pos++;
    //printf("rp =% d pos = %d\n",rp,pos);		       // PE V
    aux=0.0;
    for(int q=0;q<qmax;q++){
      aux+=w[q]*phi_r[q]*difpw[q];
      // printf(" w = %lf phi_r = %lf difpw =%lf\n",w[q],phi_r[q],difpw[q]);
    }               // PE V
    mi=indice(n0,pres,rp);
		if(std::isnan(aux*aux2) )	cout << "Float was Not a Number: ponto 3 " << aux << endl;
		
    B[mi] += aux2*aux; 			       // PE V
  }
  
  // integral 2 do caderno
	
  for(int rp=0;rp<np;rp++){
    int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
    if(rp_on_border){
      e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
      aux=0.0;
      for(int i=0;i<qmax;i++)aux+=w[i]*lambdan[0][i]*K_g_pc_n[0][i]*phi_r[i];
      mi=indice(n0,pres,rp);
      B[mi] -= aux;
    }
  }
  
  // integrais 4 e 9 do caderno
  
  for(int rp=0;rp<np;rp++){
    aux=0.0;
    for(int i=0;i<qmax;i++){
      aux+=w[i]*lambdan[iEr][i]*K_g_phi_n[iEr][pres][rp][i]*difpc[i];
    }
    mi=indice(n0,pres,rp);
    B[mi] += aux;
  }
  
  // integrais 6 e 10 do caderno
  
  for(int rp=0;rp<np;rp++){
    int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
    if(rp_on_border){
      e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
      aux=0.0;
      for(int i=0;i<qmax;i++)aux+=w[i]*difpc[i]*phi_r[i];
      mi=indice(n0,pres,rp);
      B[mi] += kln*penalty1*aux;// 07/08/08
    }
  }
  
  
  // ****************************
  // Saturation Equation - Vector
  // ****************************
  
  //integral 1
	
  for(int rs=0;rs<ns;rs++){
    int rs_on_border=e[iEr].show_ptr_stdel(sat)->is_on_border(rs,a[iEr],pos);
    if(rs_on_border){
      e[iEr].Traco_phi_1(a[iEr],sat,pos,phi_r);
      aux=0.0;
      for(int i=0;i<qmax;i++)aux+=w[i]*lambdaw[0][i]*K_g_pw_n[0][i]*phi_r[i];
      mi=indice(n0,sat,rs);
      B[mi] -= aux;
    }
  }
  
  // integrais 2 e 4
  
  for(int rs=0;rs<ns;rs++){
    aux=0.0;
    for(int i=0;i<qmax;i++){
      aux+=w[i]*lambdaw[iEr][i]*K_g_phi_n[iEr][sat][rs][i]*difpw[i];
    }
    mi=indice(n0,sat,rs);
    B[mi] += aux;
  }
  
  // integrais  3 e 5
  for(int rs=0;rs<ns;rs++){
    int rs_on_border=e[iEr].show_ptr_stdel(sat)->is_on_border(rs,a[iEr],pos);
    if(rs_on_border){
      e[iEr].Traco_phi_1(a[iEr],sat,pos,phi_r);
      
      aux=0.0;
      for(int i=0;i<qmax;i++)aux+=w[i]*difpw[i]*phi_r[i];
      mi=indice(n0,sat,rs);
      B[mi] += klw*penalty*aux;
    }
  }
  // ok 26/05/2008
  
  // ***************************************************************************
  // ********
  // Matriz *
  // ********
  // r = Pressure
  for(int rp=0;rp<np;rp++){// PE_BE
    mi=indice(n0,pres,rp);
    int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
    if(rp_on_border)e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
    
    // l = Pressure
    for(int lp=0;lp<np;lp++){ // PE_BE_PP_00
      mj=indice(n0,pres,lp);
      int lp_on_border=e[iEl].show_ptr_stdel(pres)->is_on_border(lp,a[iEl],pos);
			if(lp_on_border)e[iEl].Traco_phi_1(a[iEl],pres,pos,phi_l);
			
      aux=0.0; // inicia aux
      if(rp_on_border){
				aux1=0.0;
				for(int i=0;i<qmax;i++){
					aux1+=w[i]*lambdat[iEl][i]*K_g_phi_n[iEl][pres][lp][i]*phi_r[i];
				}
				aux-=aux1;
      }
      if(lp_on_border){
				aux1=0.0;
				for(int i=0;i<qmax;i++){
					aux1+=w[i]*lambdat[iEr][i]*K_g_phi_n[iEr][pres][rp][i]*phi_l[i];
				}
				aux += aux1;
				if(rp_on_border){
					aux1=0.0;
					for(int i=0;i<qmax;i++)aux1+=w[i]*phi_l[i]*phi_r[i];
					aux += klw*penalty*aux1;
				}
      }
			
      mx[npos(ntot,mi,mj)]+=aux;
    }
    
    // l = saturation
    for(int ls=0;ls<ns;ls++){ // PE_BE_SP_00
      mj=indice(n0,sat,ls);
      int ls_on_border=e[iEl].show_ptr_stdel(sat)->is_on_border(ls,a[iEl],pos);
			
      aux=0.0; // inicia
      
      if(ls_on_border){
				e[0].Traco_phi_1(a[iEl],sat,pos,phi_l);
				
				aux1=0.0;
				for(int q=0;q<qmax;q++){
					aux1 += w[q]*phi_l[q]*
					( d_lambdat[iEl][q]*K_g_phi_n[iEr][pres][rp][q]*difpw[q] +
					 d_lambdan[iEl][q]*K_g_phi_n[iEr][pres][rp][q]*difpc[q] +
					 lambdan[iEr][q]*d_pc[iEl][q]*K_g_phi_n[iEr][pres][rp][q] );
				}
				aux += aux1;
				
				if(rp_on_border){
					aux1=0.0;
					for(int q=0;q<qmax;q++){
						aux1 -= w[q]*phi_l[q]*phi_r[q]*
						( d_lambdat[iEl][q]*K_g_pw_n[iEl][q]
						 +d_lambdan[iEl][q]*K_g_pc_n[iEl][q]
						 +lambdan[iEl][q]*d2_pc[iEl][q]*K_g_sn_n[iEl][q]
						 -kln*penalty1*d_pc[iEl][q] );
					}
					aux += aux1;
				}
      }
      
      // integral 5 PE_BE_SP_00
      if(rp_on_border){
				aux1=0.0;
				for(int q=0;q<qmax;q++)
					aux1+=w[q]*lambdan[iEl][q]*d_pc[iEl][q]*K_g_phi_n[iEl][sat][ls][q]*phi_r[q];
				aux -= aux1;
      }
			
      mx[npos(ntot,mi,mj)]+=aux;
    }
  }
  
  // r = saturation
  for(int rs=0;rs<ns;rs++){// SE_BE
    mi=indice(n0,sat,rs);
    int rs_on_border=e[iEr].show_ptr_stdel(sat)->is_on_border(rs,a[iEr],pos);
    if(rs_on_border)e[0].Traco_phi_1(a[iEr],sat,pos,phi_r);
    
    // l = Pressure
    for(int lp=0;lp<np;lp++){ // SE_BE_PS_00
      mj=indice(n0,pres,lp);
      int lp_on_border=e[iEl].show_ptr_stdel(pres)->is_on_border(lp,a[iEl],pos);
      if(lp_on_border)e[0].Traco_phi_1(a[iEl],pres,pos,phi_l);
      aux=0.0; // inicia
      if(rs_on_border){
				//integral 1 SE_BE_PS_00
				aux1=0.0;
				for(int q=0;q<qmax;q++)
					aux1+=w[q]*lambdaw[iEl][q]*K_g_phi_n[iEl][pres][lp][q]*phi_r[q];
				aux -= aux1;
				if(lp_on_border){
					aux2=klw*penalty;
					aux1=0.0;
					for(int q=0;q<qmax;q++)
						aux1+=w[q]*phi_l[q]*(lambdaw[iEr][q]*K_g_phi_n[iEr][sat][rs][q]+aux2*phi_r[q]);
					aux += aux1;
				}
      }
			
      mx[npos(ntot,mi,mj)]+=aux;
    }
    
    for(int ls=0;ls<ns;ls++){ // SE_BE_SS_00
      mj=indice(n0,sat,ls);
      int ls_on_border=e[iEl].show_ptr_stdel(sat)->is_on_border(ls,a[iEl],pos);
      aux=0.0;// inicia
      if(ls_on_border){
				e[0].Traco_phi_1(a[iEl],sat,pos,phi_l);
				
				// integral 1 SE_BE_SS_00
				if(rs_on_border){
					aux1=0.0;
					for(int q=0;q<qmax;q++)
						aux1+=w[q]*d_lambdaw[iEl][q]*phi_l[q]*K_g_pw_n[iEl][q]*phi_r[q];
					aux -= aux1;
				}
				
				// integral 2 SE_BE_SS_00
				aux1=0.0;
				for(int q=0;q<qmax;q++)
					aux1+=w[q]*d_lambdaw[iEr][q]*phi_l[q]*K_g_phi_n[iEr][sat][rs][q]*difpw[q];
				aux += aux1;
      }
			
      mx[npos(ntot,mi,mj)]+=aux;
    }
  }
	
  // ***************************************************************************
	// ***************************************************************************
  for(int rp=0;rp<np;rp++) {
    mi=indice(n0,pres,rp);
    indx[mi] = e[iEr].map(pres,rp);
  }
  for(int rs=0;rs<ns;rs++) {
    mi=indice(n0,sat,rs);
    indx[mi] = e[iEr].map(sat,rs);
  }
	/*
	 //printf("Vetor\n");
	 for(int i=0;i<ntot;i++){
	 printf("indice  %d indx[%d] = %d  B = %lf\n",i,i,indx[i], B[i]);
	 } 
	 for(int i=0;i<iaux;i++){
	 printf("indice  %d mx = %lf \n",i,mx[i]);
	 }
	 */
	return ntot;
}
// **************************************



