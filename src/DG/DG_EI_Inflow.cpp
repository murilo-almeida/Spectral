#include "spectral.h"
#include "DG_Prob.h"

//#define indice(n0,var,i) ((n0*var)+i)
//#define npos(ntot,i,j) ((ntot*i)+j)
// ****************************************************************************
void DG_Prob::DG_EI_Inflow(const EDGE border,
                           int & count,
                           int * Ti,
                           int * Tj,
                           double * Tx,
                           double * B)
{
  //cout << "DG_EI_Inflow\n";
	// Comeco do comentario
  // Fim do comentario
	
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
//	const int n0=e[0].show_ptr_stdel(0)->nn_val();// Num de modos para a variavel 0

	int mi,mj;
	double  mx [ntot*ntot];
	//int     indx[ntot];
	
  for(int i=0;i<ntot*ntot;i++) {
    mx[i]=0.0;
  }
  
  // *******************
  //     VETOR         *
  // *******************
  // *************************
  // Pressure Equation Vector
  // *************************

  //integral 1
  pos=0;
  I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);		       // PE V
  I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);		       // PE V
  for(int ip=I0;ip<I1;ip++){					       // PE V
    int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip);	
    e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);pos++;
    
    aux=0.0;
    for(int q=0;q<qmax;q++)aux+=w[q]*phi_r[q]*lambdat[0][q]*K_g_pw_n[0][q];
    B[e[iEr].map(pres,rp)] -= aux;
  }
  
  // integrais 3 e 7 do caderno 
  
  for(int rp=0;rp<np;rp++){
    
    aux=0.0;
    for(int q=0;q<qmax;q++){
      aux+=w[q]*K_g_phi_n[iEr][pres][rp][q]*lambdat[iEr][q]*difpw[q];
    }
    
  B[e[iEr].map(pres,rp)] += aux;
  }
  
  // integrais 5 e 8 do caderno
  aux2=klt*penalty;
  I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);		       // PE V
  I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);		       // PE V
  pos=0;
  for(int ip=I0;ip<I1;ip++){					       // PE V
    int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip);		       // PE V
    e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);pos++;		       // PE V  
    aux=0.0;
    for(int q=0;q<qmax;q++)aux+=w[q]*phi_r[q]*difpw[q];                // PE V
    B[e[iEr].map(pres,rp)] += aux2*aux; 			       // PE V
  }	
  
  // integral 2 do caderno
 
  for(int rp=0;rp<np;rp++){
    int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
    if(rp_on_border){
      e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
      aux=0.0;
      for(int q=0;q<qmax;q++)aux+=w[q]*lambdan[0][q]*K_g_pc_n[0][q]*phi_r[q];
      B[e[iEr].map(pres,rp)] -= aux;
    }
  }
  
  // integrais 4 e 9 do caderno 
  
  for(int rp=0;rp<np;rp++){
    aux=0.0;
    for(int q=0;q<qmax;q++){
      aux+=w[q]*lambdan[iEr][q]*K_g_phi_n[iEr][pres][rp][q]*difpc[q];
    }
    B[e[iEr].map(pres,rp)] += aux;
  }
  
  // integrais 6 e 10 do caderno
  
  for(int rp=0;rp<np;rp++){
    int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
    if(rp_on_border){ 
      e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
      aux=0.0;
      for(int q=0;q<qmax;q++)aux+=w[q]*difpc[q]*phi_r[q];
      B[e[iEr].map(pres,rp)] += kln*penalty1*aux;// 07/08/08 
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
      for(int q=0;q<qmax;q++)aux+=w[q]*lambdaw[0][q]*K_g_pw_n[0][q]*phi_r[q];
      B[e[iEr].map(sat,rs)] -= aux;
    }
  }
  
  // integrais 2 e 4
  
  
  for(int rs=0;rs<ns;rs++){
    aux=0.0;
    for(int q=0;q<qmax;q++){
      aux+=w[q]*lambdaw[iEr][q]*K_g_phi_n[iEr][sat][rs][q]*difpw[q];
    }
    B[e[iEr].map(sat,rs)] += aux;
  }
  
  // integrais  3 e 5
  for(int rs=0;rs<ns;rs++){
    int rs_on_border=e[iEr].show_ptr_stdel(sat)->is_on_border(rs,a[iEr],pos);
    if(rs_on_border){ 
      e[iEr].Traco_phi_1(a[iEr],sat,pos,phi_r);
      
      aux=0.0;
      for(int q=0;q<qmax;q++)aux+=w[q]*difpw[q]*phi_r[q];
      B[e[iEr].map(sat,rs)] += klw*penalty*aux;
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
	for(int q=0;q<qmax;q++){
	  aux1+=w[q]*lambdat[iEl][q]*K_g_phi_n[iEl][pres][lp][q]*phi_r[q];
	}
	aux-=aux1;
      }
      if(lp_on_border){
	aux1=0.0;
	for(int q=0;q<qmax;q++){
	  aux1+=w[q]*lambdat[iEr][q]*K_g_phi_n[iEr][pres][rp][q]*phi_l[q];
	}
	aux += aux1;
	if(rp_on_border){
	  aux1=0.0;
	  for(int q=0;q<qmax;q++)aux1+=w[q]*phi_l[q]*phi_r[q];
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
    
    for(int lp=0;lp<np;lp++) {
      mj=indice(n0,pres,lp);
      Ti[count]=e[iEr].map(pres,rp);
      Tj[count]=e[iEl].map(pres,lp);
      Tx[count]= mx[npos(ntot,mi,mj)];
      count++;
    }

    for(int ls=0;ls<ns;ls++) {
      mj=indice(n0,sat,ls);
      Ti[count]=e[iEr].map(pres,rp);
      Tj[count]=e[iEl].map(sat,ls);
      Tx[count]= mx[npos(ntot,mi,mj)];
      count++;
    }
  }
  for(int rs=0;rs<ns;rs++) {
    mi=indice(n0,sat,rs);

    for(int lp=0;lp<np;lp++) {
      mj=indice(n0,pres,lp); 
      Ti[count]=e[iEr].map(sat,rs);
      Tj[count]=e[iEl].map(pres,lp);
      Tx[count]= mx[npos(ntot,mi,mj)];
      count++;
    }
  
    for(int ls=0;ls<ns;ls++) {
      mj=indice(n0,sat,ls);   
      Ti[count]=e[iEr].map(sat,rs);
      Tj[count]=e[iEl].map(sat,ls);
      Tx[count]= mx[npos(ntot,mi,mj)];
      count++;
    }
  }
};
// alterado em 21/09/2011
// alterado em 13/09/2012
// alterado em 30/11/2013
// alterado em 23/05/2015 Alguns loops i<qmax foram convertidos para q<qmax
