//
//  DG_EI_Interior.cpp
//  SDG
//
//  Created by Murilo Almeida on 11/30/13.
//  Copyright (c) 2013 Murilo Almeida. All rights reserved.
//

#include "spectral.h"
#include "DG_Prob.h"
int DG_Prob::DG_EI_Interior(const EDGE border,
														double * mx,
														double * B,
														int    * indx)
{
	// *************************************************************************
#include "DG_EI_Header.h" // Inicia variaveis locais comuns a DG_EI
	// *************************************************************************
	
	
	// **************************************************************************
	// Interior Edges
	// **************************************************************************
	
	double difpw[qmax],difpc[qmax];
	for (int q=0;q<qmax;q++) {
		difpw[q]=pw[0][q]-pw[1][q];
		difpc[q]=pc[0][q]-pc[1][q];
	}
	
	// *************************************************************************
	
	int ns[2],np[2],nnparc[4];
	int ntot;
	ns[0]=e[0].show_ptr_stdel(sat)->nn_val();
	np[0]=e[0].show_ptr_stdel(pres)->nn_val();
	ns[1]=e[1].show_ptr_stdel(sat)->nn_val();
	np[1]=e[1].show_ptr_stdel(pres)->nn_val();
	
    
	nnparc[0]=0;
  nnparc[1]=          e[0].show_ptr_stdel(0)->nn_val();
  nnparc[2]=nnparc[1]+e[0].show_ptr_stdel(1)->nn_val();
  nnparc[3]=nnparc[2]+e[1].show_ptr_stdel(0)->nn_val();
  ntot     =nnparc[3]+e[1].show_ptr_stdel(1)->nn_val();
	/*
	nnparc[0]=0;
	nnparc[1]=ns[0];
	nnparc[2]=nnparc[1]+np[0];
	nnparc[3]=nnparc[2]+ns[1];
	ntot     =nnparc[3]+np[1];
	*/
    
	int mi,mj;
	int iaux=ntot*ntot;
	
	for(int i=0; i<iaux; ++i) {
		mx[i]=0.0;
	}
	for(int i=0; i<ntot; ++i) {
		B[i]=0.0;
	}
	
	// **************************************************************************
	// Pressure Equation Vector
	// **************************************************************************
	double aux3,aux4,aux5,aux6;
	double res0[qmax],res1[qmax],res2[qmax];
	
	int temp;							        // PE V
	// PE V
	//integral (1)						        // PE V
	
	for(int i=0; i<qmax; i++)
		res0[i]=lambdat[0][i]*K_g_pw_n[0][i]+lambdat[1][i]*K_g_pw_n[1][i];
	// dobro de {lambdat K nabla pw . ne}
	temp=1;							        // PE V
	for(iEr=0;iEr<2;iEr++) {					        // PE V
		pos=0;							        // PE V
		int I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);	        // PE V
		int I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);	        // PE V
		for(int ip=I0;ip<I1;ip++) {					        // PE V
			int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip);		        // PE V
			
			e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
			pos++;				       // PE V
			aux=0.0;
			for(int q=0;q<qmax;q++) aux+=w[q]*res0[q]*phi_r[q];// PE V
      //B[e[iEr].map(pres,rp)] -= 0.5*aux*temp;
      mi=indice_interior(nnparc,iEr,pres,rp);
      B[mi] -= 0.5*aux*temp;		        // PE V
		}								        // PE V
		temp*=-1;							        // PE V
	}								        // PE V
	// PE V
	// integral (2)						        // PE V
	// PE V
	for(int i=0; i<qmax; i++){
		res0[i]=0.0;
		for(int ndir=0;ndir<ndim;ndir++) {
			res0[i]+=lambdan[ndir][i]*K_g_pc_n[ndir][i];
		}
	}
	temp=1;							         // PE V
	for(iEr=0;iEr<2;iEr++) {					         // PE V
		pos=0; 						                 // PE V
		int I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);		 // PE V
		int I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);		 // PE V
		for(int ip=I0;ip<I1;ip++) {						 // PE V
			int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip); 		 // PE V
			// PE V
			e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
			pos++;				                                 // PE V
			aux=0.0;
			for(int q=0;q<qmax;q++) aux+=w[q]*res0[q]*phi_r[q];
      // PE V
      //B[e[iEr].map(pres,rp)] -= 0.5*aux*temp;
      mi=indice_interior(nnparc,iEr,pres,rp);
      B[mi] -= 0.5*aux*temp;				 // PE V
		}									 // PE V
		temp*=-1;								 // PE V
	}									 // PE V
	// PE V
	// integrais (3) e (4)						 // PE V
	
	for(iEr=0;iEr<2;iEr++) {						 // PE V
		for(int rp=0;rp<np[iEr];rp++) {					 // PE V
			aux=0.0;
			for(int q=0;q<qmax;q++)
				aux+=w[q]*K_g_phi_n[iEr][pres][rp][q]*(lambdat[iEr][q]*difpw[q]+
																							 lambdan[iEr][q]*difpc[q]);// PE V
      //B[e[iEr].map(pres,rp)] += 0.5*aux;
      mi=indice_interior(nnparc,iEr,pres,rp);
      B[mi] += 0.5*aux;				 // PE V
		}									 // PE V
	}									 // PE V
	// PE V
	// integral (5)							 // PE V
	
	temp=1;								 // PE V
	for(iEr=0;iEr<2;iEr++) {
		pos=0;						                 // PE V
		int I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);		 // PE V
		int I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);		 // PE V
		for(int ip=I0;ip<I1;ip++) {						 // PE V
			int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip);			 // PE V
			e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
			pos++;				                                 // PE V
			aux=0.0;
			for(int q=0;q<qmax;q++)aux+=w[q]*phi_r[q]*difpw[q]; 		 // PE V
      //B[e[iEr].map(pres,rp)] += klt*penalty*aux*temp;
      mi=indice_interior(nnparc,iEr,pres,rp);
      B[mi] += klt*penalty*aux*temp; 			 // PE V
		}									 // PE V
		temp*=-1;								 // PE V
	}									 // PE V
	// PE V
	// integral (6)							 // PE V
	// PE V
	temp=1;								 // PE V
	for(iEr=0;iEr<2;iEr++) {						 // PE V
		pos=0;								 // PE V
		int I0=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]);		 // PE V
		int I1=(e[iEr].show_ptr_stdel(pres))->show_bmapi(a[iEr]+1);		 // PE V
		for(int ip=I0;ip<I1;ip++) {						 // PE V
			int rp=(e[iEr].show_ptr_stdel(pres))->show_bmapv(ip);			 // PE V
			// PE V
			e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
			pos++;				                                 // PE V
			aux=0.0;
			for(int q=0;q<qmax;q++)aux+=w[q]*phi_r[q]*difpc[q];		 // PE V
      // B[e[iEr].map(pres,rp)] += kln*penalty1*aux*temp;// 07/08/08 	 // PE V
      mi=indice_interior(nnparc,iEr,pres,rp);
      B[mi]  += kln*penalty1*aux*temp;
		}									 // PE V
		temp*=-1;								 // PE V
	}                                                                      // PE V
	
	// ****************************                                        // SE V
	// Saturation Equation Vector *					 // SE V
	// ****************************					 // SE V
	// SE V
	// integral (1)							 // SE V
	// SE V
	for(int q=0;q<qmax;q++){
		res2[q]=0.0;
		for(int ndir=0;ndir<ndim;ndir++){
			res2[q]+=lambdaw[ndir][q]*K_g_pw_n[ndir][q];
		}
	}
	temp=1;								 // SE V
	for(iEr=0;iEr<2;iEr++) {						 // SE V
		pos=0;								 // SE V
		int I0=(e[iEr].show_ptr_stdel(sat))->show_bmapi(a[iEr]);		 // SE V
		int I1=(e[iEr].show_ptr_stdel(sat))->show_bmapi(a[iEr]+1);		 // SE V
		for(int ip=I0;ip<I1;ip++) {						 // SE V
			int rs=(e[iEr].show_ptr_stdel(sat))->show_bmapv(ip); 			 // SE V
			// SE V
			e[iEr].Traco_phi_1(a[iEr],sat,pos,phi_r);
			pos++;				                                 // SE V
			aux=0.0;
			for(int q=0;q<qmax;q++)aux+=w[q]*phi_r[q]*res2[q];
      //B[e[iEr].map(sat,rs)] -= 0.5*aux*temp;				 // SE V
      mi=indice_interior(nnparc,iEr,sat,rs);
      B[mi]  -= 0.5*aux*temp;
		}									 // SE V
		temp*=-1;								 // SE V
	}									 // SE V
	// SE V
	// integral (2)							 // SE V
	// SE V
	for(iEr=0;iEr<2;iEr++) {						 // SE V
		for(int rs=0;rs<ns[iEr];rs++) {					 // SE V
			aux=0.0;
			for(int q=0;q<qmax;q++)
				aux+=w[q]*lambdaw[iEr][q]*K_g_phi_n[iEr][sat][rs][q]*difpw[q];	 // SE V
      //B[e[iEr].map(sat,rs)] += 0.5*aux;
      mi=indice_interior(nnparc,iEr,sat,rs);
      B[mi]  += 0.5*aux;				 // SE V
		}									 // SE V
	}									 // SE V
	// SE V
	// integral (3)							 // SE V
	// SE V
	temp=1;								 // SE V
	for(iEr=0;iEr<2;iEr++) {						 // SE V
		pos=0;								 // SE V
		int I0=(e[iEr].show_ptr_stdel(sat))->show_bmapi(a[iEr]);		 // SE V
		int I1=(e[iEr].show_ptr_stdel(sat))->show_bmapi(a[iEr]+1);		 // SE V
		for(int ip=I0;ip<I1;ip++) {						 // SE V
			int rs=(e[iEr].show_ptr_stdel(sat))->show_bmapv(ip); 			 // SE V
			// SE V
			e[iEr].Traco_phi_1(a[iEr],sat,pos,phi_r);pos++;			 // SE V
			aux=0.0;
			for(int q=0;q<qmax;q++)
				aux+=w[q]*phi_r[q]*difpw[q];
      // B[e[iEr].map(sat,rs)] += klw*penalty*aux*temp;
      mi=indice_interior(nnparc,iEr,sat,rs);
      B[mi]  += klw*penalty*aux*temp;			 // SE V
		}									 // SE V
		temp*=-1;								 // SE V
	}                                                                      // SE V
	
	// **************************************************************************
	// ********
	// Matriz *
	// ********
	// **************************************************************************
	// Pressure equation
	// **************************************************************************
	
	// ***************************************************************************
	for(iEr=0;iEr<2;iEr++) {// iEr                                      // PE_IE_P
		int sr=Sinal(iEr);						      // PE_IE_P
		// PE_IE_P
		// Pressure							      // PE_IE_P
		// PE_IE_P
		for(int rp=0;rp<np[iEr];rp++) {// PE_IE  rp			      // PE_IE_P
			
			int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
			if(rp_on_border)
				e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
      
      for(iEl=0;iEl<2;iEl++) {//  iEl				      // PE_IE_P
				int sl=Sinal(iEl);					      // PE_IE_P
				// PE_IE_P
				// Pressure						      // PE_IE_P
				for(int lp=0;lp<np[iEl];lp++) { // PE_IE_P lp		      // PE_IE_P
					
					int lp_on_border=e[iEl].show_ptr_stdel(pres)->is_on_border(lp,a[iEl],pos);
					if(lp_on_border)
						e[iEl].Traco_phi_1(a[iEl],pres,pos,phi_l);
					
					// integral 1						      // PE_IE_P
					// PE_IE_P
					if(rp_on_border) {					      // PE_IE_P
						aux1=0.0;
						for(int q=0;q<qmax;q++)
							aux1+=w[q]*lambdat[iEl][q]*K_g_phi_n[iEl][pres][lp][q]*phi_r[q];
					}							      // PE_IE_P
					else aux1=0.0;					      // PE_IE_P
					// PE_IE_P
					// integral 3						      // PE_IE_P
					// PE_IE_P
					if(lp_on_border) {					      // PE_IE_P
						aux2=0.0;
						for(int q=0;q<qmax;q++)
							aux2+=w[q]*lambdat[iEr][q]*K_g_phi_n[iEr][pres][rp][q]*phi_l[q];
						// integral 5					      // PE_IE_P
						if(rp_on_border) {				              // PE_IE_P
							aux3=0.0;
							for(int q=0;q<qmax;q++)aux3+=w[q]*phi_l[q]*phi_r[q];    // PE_IE_P
						}							      // PE_IE_P
						else aux3=0.0;					      // PE_IE_P
					}							      // PE_IE_P
					else {						      // PE_IE_P
						aux2=0.0;						      // PE_IE_P
						aux3=0.0;						      // PE_IE_P
					}							      // PE_IE_P
					// PE_IE_P
					aux = -0.5*sr*aux1+0.5*sl*aux2+sl*sr*klt*penalty*aux3;      // PE_IE_P
					{
						mi=indice_interior(nnparc,iEr,pres,rp);
						mj=indice_interior(nnparc,iEl,pres,lp);
						mx[npos(ntot,mi,mj)]+=aux;			              // PE_IE_P
					}						              // PE_IE_P
				}                                                             // PE_IE_P
				
				// Saturation
				
				for(int ls=0;ls<ns[iEl];ls++) { // PE_IE_S ls                 // PE_IE_S
					// PE_IE_S
					int ls_on_border=e[iEl].show_ptr_stdel(sat)->is_on_border(ls,a[iEl],pos);
					if(ls_on_border)
						e[iEl].Traco_phi_1(a[iEl],sat,pos,phi_l);                 // PE_IE_S
					aux=0.0;						      // PE_IE_S
					// PE_IE_S
					if(rp_on_border) {					      // PE_IE_S
						// PE_IE_S
						if(ls_on_border) {					      // PE_IE_S
							// PE_IE_S
							// integrais 1,2,3 				      // PE_IE_S
							aux1=0.0;
							for(int q=0;q<qmax;q++) {
								aux1+=w[q]*phi_l[q]*phi_r[q]*
								(d_lambdat[iEl][q]*K_g_pw_n[iEl][q]-
								 d_lambdan[iEl][q]*K_g_pc_n[iEl][q]-
								 lambdan[iEl][q]*d2_pc[iEl][q]*K_g_sn_n[iEl][q]);
							}
							aux -= 0.5*sr*aux1;				      // PE_IE_S
						}							      // PE_IE_S
						// PE_IE_S
						// integral 4					      // PE_IE_S
						aux2=0.0;
						for(int q=0;q<qmax;q++) {
							aux2+=w[q]*lambdan[iEl][q]*d_pc[iEl][q]*K_g_phi_n[iEl][sat][ls][q]*phi_r[q];
						}
						aux -= 0.5*sr*aux2;					      // PE_IE_S
					}							      // PE_IE_S
					// PE_IE_S
					if(ls_on_border) {					      // PE_IE_S
						// PE_IE_S
						if(iEl==iEr) {					      // PE_IE_S
							// integral (7)					      // PE_IE_S
							aux3=0.0;
							for(int q=0;q<qmax;q++)
								aux3+=w[q]*d_lambdat[iEl][q]*phi_l[q]*
								K_g_phi_n[iEr][pres][rp][q]*difpw[q];		      // PE_IE_S
							aux += 0.5*aux3;					      // PE_IE_S
							// PE_IE_S
							// integral (9)					      // PE_IE_S
							aux4=0.0;
							for(int q=0;q<qmax;q++)
								aux4+=w[q]*d_lambdan[iEl][q]*phi_l[q]*K_g_phi_n[iEr][pres][rp][q]*difpc[q];
							aux += 0.5*aux4;					      // PE_IE_S
						}// Integrais que so existem quando iEl=iEr		      // PE_IE_S
						// PE_IE_S
						// integral (5)					      // PE_IE_S
						
						aux5=0.0;
						for(int q=0;q<qmax;q++)
							aux5+=w[q]*lambdan[iEr][q]*K_g_phi_n[iEr][pres][rp][q]*
							d_pc[iEl][q]*phi_l[q];
						aux += 0.5*aux5*sl;					      // PE_IE_S
						// PE_IE_S
						// integral (6)					      // PE_IE_S
						if(rp_on_border) {					      // PE_IE_S
							aux6=0.0;
							for(int q=0;q<qmax;q++)
								aux6+=w[q]*d_pc[iEl][q]*phi_l[q]*phi_r[q];	      // PE_IE_S
							aux += kln*penalty1*sl*sr*aux6;// 07/08/08	      // PE_IE_S
						}							      // PE_IE_S
					}							      // PE_IE_S
					/*if(aux!=0.0)*/
					{					                      // PE_IE_S
						mi=indice_interior(nnparc,iEr,pres,rp);
						mj=indice_interior(nnparc,iEl,sat,ls);
						mx[npos(ntot,mi,mj)]+=aux;						      // PE_IE_S
					}							      // PE_IE_S
				}							      // PE_IE_S
      }								      // PE_IE_S
		}								      // PE_IE_S
	}                                                                   // PE_IE_S
	
	// **************************************************************************
	// Saturation Equation
	// **************************************************************************
	// SE_IE_P
	for(iEr=0;iEr<2;iEr++) {					      // SE_IE_P
		int sr=Sinal(iEr);						      // SE_IE_P
		// SE_IE_P
		// Saturation						      // SE_IE_P
		// SE_IE_P
		for(int rs=0;rs<ns[iEr];rs++) {// SE_IE  rs			      // SE_IE_P
			// SE_IE_P
			int rs_on_border=e[iEr].show_ptr_stdel(sat)->is_on_border(rs,a[iEr],pos);
			if(rs_on_border)
				e[iEr].Traco_phi_1(a[iEr],sat,pos,phi_r);		      // SE_IE_P
      for(iEl=0;iEl<2;iEl++) {// <==================================  // SE_IE_P
				int sl=Sinal(iEr);					      // SE_IE_P
				// SE_IE_P
				// Derivative of Pressure				      // SE_IE_P
				// SE_IE_PS 						      // SE_IE_P
				// SE_IE_P
				for(int lp=0;lp<np[iEl];lp++) {				      // SE_IE_P
					// SE_IE_P
					int lp_on_border=e[iEl].show_ptr_stdel(pres)->is_on_border(lp,a[iEl],pos);
					if(lp_on_border)
						e[iEl].Traco_phi_1(a[iEl],pres,pos,phi_l);		      // SE_IE_P
					aux=0.0;						      // SE_IE_P
					// SE_IE_P
					if (rs_on_border) {					      // SE_IE_P
						//integral 1					      // SE_IE_P
						aux1=0.0;
						for(int q=0;q<qmax;q++)
							aux1+=w[q]*lambdaw[iEl][q]*K_g_phi_n[iEl][pres][lp][q]*phi_r[q];
						aux -= 0.5*sr*aux1;					      // SE_IE_P
					}							      // SE_IE_P
					// SE_IE_P
					if (lp_on_border) {					      // SE_IE_P
						// SE_IE_P
						// integral 2 					      // SE_IE_P
						aux1=0.0;
						for(int q=0;q<qmax;q++)
							aux1+=w[q]*lambdaw[iEr][q]*K_g_phi_n[iEr][sat][rs][q]*phi_l[q];
						aux += 0.5*sl*aux1;					      // SE_IE_P
						// SE_IE_P
						if (rs_on_border) {					      // SE_IE_P
							// SE_IE_P
							// integral 3					      // SE_IE_P
							aux1=0.0;
							for(int q=0;q<qmax;q++)
								aux1+=w[q]*phi_l[q]*phi_r[q];
							aux += klw*penalty*sl*sr*aux1;			      // SE_IE_P
						}							      // SE_IE_P
					}							      // SE_IE_P
					/*if(aux!=0.0)*/
					{
						mi=indice_interior(nnparc,iEr,sat,rs);
						mj=indice_interior(nnparc,iEl,pres,lp);
						mx[npos(ntot,mi,mj)]+=aux;						      // SE_IE_P
					}							      // SE_IE_P
				}                                                             // SE_IE_P
				
				// Derivative of Saturation
				
				// SE_IE_SS_00                                                // SE_IE_S
				pos=0; 							      // SE_IE_S
				int I0=(e[iEl].show_ptr_stdel(sat))->show_bmapi(a[iEl]);	      // SE_IE_S
				int I1=(e[iEl].show_ptr_stdel(sat))->show_bmapi(a[iEl]+1);	      // SE_IE_S
				for(int ip=I0;ip<I1;ip++) {				      // SE_IE_S
					int ls=(e[iEl].show_ptr_stdel(sat))->show_bmapv(ip); 	      // SE_IE_S
					// SE_IE_S
					e[iEl].Traco_phi_1(a[iEl],sat,pos,phi_l);
					pos++;		      // SE_IE_S
					
					for(int qtemp=0.0;qtemp<qmax;qtemp++) {
						res1[qtemp]=0.0;
						for(int ndir=0;ndir<ndim;ndir++){
							res1[qtemp]+=Kgpw[iEl][ndir][qtemp]*n_e[ndir];
						}
					}
					//prodgn(qmax,Kgpw[iEl][0],Kgpw[iEl][1],n_e,res1);	      // SE_IE_S
					// SE_IE_S
					aux=0.0;	  					      // SE_IE_S
					// SE_IE_S
					if (rs_on_border) {					      // SE_IE_S
						// SE_IE_S
						//integral 1 					      // SE_IE_S
						// SE_IE_S
						aux1=0.0;
						for(int q=0;q<qmax;q++)
							aux1+=w[q]*d_lambdaw[iEl][q]*K_g_pw_n[iEl][q]*phi_l[q]*phi_r[q];
						aux -= 0.5*sr*aux1;					      // SE_IE_S
					}							      // SE_IE_S
					// SE_IE_S
					if(iEl==iEr) {					      // SE_IE_S
						// integral 2					      // SE_IE_S
						aux1=0.0;
						for(int q=0;q<qmax;q++)
							aux1+=w[q]*d_lambdaw[iEl][q]*phi_l[q]*
							K_g_phi_n[iEr][sat][rs][q]*difpw[q];
						aux += 0.5*aux1;					      // SE_IE_S
					}							      // SE_IE_S
					/*if(aux!=0.0)*/
					{					                      // SE_IE_S
						mi=indice_interior(nnparc,iEr,sat,rs);
						mj=indice_interior(nnparc,iEl,sat,ls);
						mx[npos(ntot,mi,mj)]+=aux;						      // SE_IE_S
					}							      // SE_IE_S
				}							      // SE_IE_S
      }								      // SE_IE_S
		}								      // SE_IE_S
	}
	
	//printf(" Ate Aqui\n");                                                              // SE_IE_S
	// Ate aqui esta ok 22/05/2008
	// ******************************************
	// ******************************************
	for(iEr=0;iEr<2;iEr++) {
		for (int rp=0;rp<np[iEr];rp++) {
			mi = indice_interior(nnparc,iEr,pres,rp);
			indx[mi] = e[iEr].map(pres,rp);
		}
		for (int rs=0;rs<ns[iEr];rs++) {
			mi = indice_interior(nnparc,iEr,sat,rs);
			indx[mi] = e[iEr].map(sat,rs);
		}
	}
	return ntot;
}
