#include "spectral.h"
#include "DG_Prob.h"
// ****************************************************************************
// Numeracao das integrais de acordo com o caderno DG_IG
// ****************************************************************************

void DG_Prob::DG_EI_IG(const EDGE border,
                       const double sigma,
                       const double sigma1,
                       const double beta, int & count,
                       int * Ti,int * Tj, double * Tx,
                       double * B,
                       const int qmax,
                       const double w0[])
{
 
  double aux,aux1,aux2,aux3;
  int m,lado;
  MyElem e[2];
 
  int tipo=border.tipo;

  e[0]=el[border.elemento[0]];
  const int ndim=e[0].show_ptr_stdel(0)->ndim_val();// numero de dimensooes

  double w[qmax];
  //Gauss_Jacobi_parameters(qmax,0.0,0.0,Jac,w,Dtemp);
  double length;
  if(ndim==1){ 
    length=1.0;
  }
  else {
    // Gauss_Jacobi_parameters(qmax,0.0,0.0,Jac,w,Dtemp);
    length=border.comprimento; 
  }
  // *******************************************************
  // ******************************************************* 
  for(int q=0;q<qmax;q++)w[q]=w0[q]*length/2.0;// salva em w = w*length/2.0
  // *******************************************************
  const double penalty =sigma/pow(length,beta);
  const double penalty1 =sigma1/pow(length,beta);
  //printf("sigma = %11.4e beta = %11.4e length = %11.4e penalty = %11.4e\n",sigma, beta,length,penalty);

  double n_e[2];

  double perm[2][3];
  int a[2];
  double sn[2][qmax],pw[2][qmax];
  double lambdaw[2][qmax],lambdan[2][qmax],lambdat[2][qmax];
  double d_lambdaw[2][qmax],d_lambdan[2][qmax],d_lambdat[2][qmax];
  double pc[2][qmax],d_pc[2][qmax],d2_pc[2][qmax];
  double phi_r[qmax];
  double phi_l[qmax];

  double mun=fluids.show_mun();
  double muw=fluids.show_muw();
 
  double res1[qmax],res2[qmax],res3[qmax],res0[qmax];
  int ns[2],np[2];
  int iEr, iEl;
  double kn,kln,klt; 

  for(int i=0; i<ndim;i++){
    n_e[i]=border.normal[i];
  }
  
  lado=0;
  e[0].perm_val(perm[0]);
  a[0]=border.num_local[0];

  e[0].Traco_sn(a[0],sn[0]);
  e[0].Traco_pw(a[0],pw[0]);
 
  for(m=0;m<qmax;m++){
    // Calculo dos escalares
    double aux=sn[lado][m];
    lambdan[lado][m]=fluids.Krn(aux)/mun;
    lambdaw[lado][m]=fluids.Krw(aux)/muw;
    lambdat[lado][m]=lambdan[lado][m]+lambdaw[lado][m];
    d_lambdan[lado][m]=fluids.dkrn(aux)/mun;
    d_lambdaw[lado][m]=fluids.dkrw(aux)/muw;
    d_lambdat[lado][m]=d_lambdan[lado][m]+d_lambdaw[lado][m];
    pc[lado][m]=fluids.pressao_capilar(aux);
    d_pc[lado][m]=fluids.dpc(aux);
    d2_pc[lado][m]=fluids.d2pc(aux);
  }
 
  aux2=0.0;
  for(int i=0;i<ndim;i++){
    aux=perm[0][i]*n_e[i];
    aux2+=aux*aux;
  }
  kn=sqrt(aux2); // K pre-fator
 
  //prodev(qmax,kn,lambdan[0],res0);
  //integral(qmax,w,res0,aux);
  //kln=aux/length;

  //prodev(qmax,kn,lambdat[0],res0);
  //integral(qmax,w,res0,aux);
  //klt=aux/length;
 
  kln=1.0; //kn/mun;
  klt=1.0; //kln+kn/muw;

  e[lado].Traco_Kgrad_sn(a[lado],Kgsn[lado]);
  e[lado].Traco_Kgrad_pw(a[lado],Kgpw[lado]);
  e[lado].Traco_Kgrad_pc(a[lado],Kgpc[lado]);
  
  np[0]=e[0].show_ptr_stdel(pres)->nn_val();
  ns[0]=e[0].show_ptr_stdel(sat)->nn_val();
  
  // Echo dos dados
  
 //   fprintf(fout1,"%4d %4d %4d ",index,border.Na,border.Nb);index++;
 //   int e0=border.elemento[0];
 //   int s0=border.sinal[0];
 //   fprintf(fout1,"   %4d %4d %5d ",e0,a[0],s0);
 //   if(tipo==2){
 //     e0=border.elemento[1];
 //     s0=border.sinal[1];
 //     a[1]=border.num_local[1];
 //     fprintf(fout1,"   %4d %4d %5d  ",e0,a[1],s0);
 //   }
 //   else fprintf(fout1,"          EXTERIOR  ");
 //   
 //   fprintf(fout1," %11.4e %11.4e\n",n_e[0],n_e[1]);
 //   printf("Inside EdgeIntegrals: tipo = %2d \n",tipo);
  // Fim do Echo dos dados
  
  // **************************************************************************
  // Interior Edges
  // **************************************************************************

  if(tipo==2){ // tipo==2  Interior Edges
    lado=1;
    a[lado]=border.num_local[lado];
    e[lado]=el[border.elemento[lado]];
    e[lado].perm_val(perm[lado]);
    e[lado].Traco_sn(a[lado],sn[lado]); 
    e[lado].Traco_pw(a[lado],pw[lado]);

    for(m=0;m<qmax;m++){
      // Calculo dos escalares
      double aux=sn[lado][m];
      lambdan[lado][m]=fluids.Krn(aux)/mun;
      lambdaw[lado][m]=fluids.Krw(aux)/muw;
      lambdat[lado][m]=lambdan[lado][m]+lambdaw[lado][m];
      d_lambdan[lado][m]=fluids.dkrn(aux)/mun;
      d_lambdaw[lado][m]=fluids.dkrw(aux)/muw;
      d_lambdat[lado][m]=d_lambdan[lado][m]+d_lambdaw[lado][m];
      pc[lado][m]=fluids.pressao_capilar(aux);
      d_pc[lado][m]=fluids.dpc(aux);
      d2_pc[lado][m]=fluids.d2pc(aux);
    }
    
    aux2=0.0;
    for(int i=0;i<ndim;i++){
      aux=perm[1][i]*n_e[i];
      aux2+=aux*aux;
    }
    kn=sqrt(aux2); // K pre-fator
    if(std::isnan(aux2))cout << "DG_EI_IG  Float was Not a Number: aux2 " << aux2 << endl;
    //prodev(qmax,kn,lambdan[1],res0);
    //integral(qmax,w,res0,aux);
    /*    
	  aux=kn/mun;
    kln = kln < aux ? kln : aux;
     if(std::isnan(kln))cout << "DG_EI_IG  Float was Not a Number: kln " << kln << endl;
    //prodev(qmax,kn,lambdat[1],res0);
    //integral(qmax,w,res0,aux);
    aux=kln+kn/muw;
    klt = klt < aux ? klt : aux;
    */
    kln=1.0;
    klt=1.0; 

    e[lado].Traco_Kgrad_sn(a[lado],Kgsn[lado]);
    e[lado].Traco_Kgrad_pw(a[lado],Kgpw[lado]);
    e[lado].Traco_Kgrad_pc(a[lado],Kgpc[lado]);

    np[1]=e[1].show_ptr_stdel(pres)->nn_val();
    ns[1]=e[1].show_ptr_stdel(sat)->nn_val();
   
    // ************************************************************************
    // Pressure Equation Vector
    // ************************************************************************
    int temp;

    // integral 2
    
    for(int qtemp=0;qtemp<qmax;qtemp++) {
      res0[qtemp]=0.0;
      for(int ndir=0;ndir<ndim;ndir++){
	res0[qtemp]+=Kgpc[0][ndir][qtemp]*n_e[ndir];
      }
    }
    
    //prodgn(qmax,Kgpc[0][0],Kgpc[0][1],n_e,res0);
    
    for(int qtemp=0;qtemp<qmax;qtemp++) {
      res1[qtemp]=0.0;
      for(int ndir=0;ndir<ndim;ndir++){
	res1[qtemp]+=Kgpc[1][ndir][qtemp]*n_e[ndir];
      }
    }
    
    //prodgn(qmax,Kgpc[1][0],Kgpc[1][1],n_e,res1);
    
    prodvv(qmax,lambdan[0],res0,res0);
    prodvv(qmax,lambdan[1],res1,res1);
    
    //somavv(qmax,res0,res1,res0);
    for(int i=0;i<qmax;i++)res0[i]+=res1[i];
    temp=1;
    for(iEr=0;iEr<2;iEr++){
      for(int rp=0;rp<np[iEr];rp++){
//	int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
	//if(rp_on_border){
	//e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
	e[iEr].Traco_phi(a[iEr],pres,rp,phi_r);
	prodvv(qmax,res0,phi_r,res1);
	integral(qmax,w,res1,aux);
	B[e[iEr].map(pres,rp)] += 0.5*aux*temp;
	//}
      }
      temp*=-1;
    }
    // OK 22/05/2008
    // integral 4

    subtvv(qmax,pc[0],pc[1],res0);// salva res0
    for(iEr=0;iEr<2;iEr++){
      for(int rp=0;rp<np[iEr];rp++){
        //e[iEr].Traco_grad_phi(a[iEr],pres,rp,gphi_r);
	
        e[iEr].Traco_Kgrad_phi_n(a[iEr],pres,rp,res2);

        prod3v(qmax,lambdan[iEr],res2,res0,res1);
        integral(qmax,w,res1,aux);
        if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 3 " << aux << endl;
        B[e[iEr].map(pres,rp)] -= 0.5*aux;
      }
    }
    // OK 22/05/2008
    // integral 6

    temp=1;
    for(iEr=0;iEr<2;iEr++){
      for(int rp=0;rp<np[iEr];rp++){
	// int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
	//if(rp_on_border){
	// e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
	e[iEr].Traco_phi(a[iEr],pres,rp,phi_r);
	prodvv(qmax,res0,phi_r,res1);
	integral(qmax,w,res1,aux);
	if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 4 " << aux << endl;
	B[e[iEr].map(pres,rp)] -= kln*penalty1*aux*temp;// 07/08/08
	//}
      }
      temp*=-1;
    }
    // OK 22/05/2008
    
    // ********
    // Matriz *
    // ********
    // ************************************************************************
    // Pressure equation
    // ************************************************************************
    
    for(iEr=0;iEr<2;iEr++){// iEr
      int sr=Sinal(iEr);
      
      // Pressure
      
      for(int rp=0;rp<np[iEr];rp++){// PE_IE  rp
	//int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
	//if(rp_on_border){
	//e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
	e[iEr].Traco_phi(a[iEr],pres,rp,phi_r);
	//e[iEr].Traco_grad_phi(a[iEr],pres,rp,gphi_r);
	e[iEr].Traco_Kgrad_phi_n(a[iEr],pres,rp,res0);
	
	for(iEl=0;iEl<2;iEl++){//  iEl
	  int sl=Sinal(iEl);
	  
	  // Pressure
	  for(int lp=0;lp<np[iEl];lp++){ // PE_IE_P lp
//	    int lp_on_border=e[iEl].show_ptr_stdel(pres)->is_on_border(lp,a[iEl],pos);
	    //if(lp_on_border){
	    //e[iEl].Traco_phi_1(a[iEl],pres,pos,phi_l);
	    e[iEl].Traco_phi(a[iEl],pres,lp,phi_l);
	   // e[iEl].Traco_grad_phi(a[iEl],pres,lp,gphi_l);
	    e[iEl].Traco_Kgrad_phi_n(a[iEl],pres,lp,res1);
	    
	    // integral 1
	    
	    prod3v(qmax,lambdat[iEl],res1,phi_r,res1);
	    integral(qmax,w,res1,aux1);
	    
	    // integral 3
	    
	    prod3v(qmax,lambdat[iEr],res0,phi_l,res2);
	    integral(qmax,w,res2,aux2);
	    
	    // integral 5
	    
	    prodvv(qmax,phi_l,phi_r,res3);
	    integral(qmax,w,res3,aux3);
	    
	    aux = -0.5*sr*aux1+0.5*sl*aux2+sl*sr*klt*penalty*aux3;
	    if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 4 " << aux << endl;
	    Ti[count]=e[iEr].map(pres,rp);
	    Tj[count]=e[iEl].map(pres,lp);
	    Tx[count]=aux;
	    count++;
	    //}
	  }
	  // OK 22/05/2008  
	}
	//}
      }
    }
  }
  
  // Ate aqui esta ok 22/05/2008  
  
  // FIM do if(tipo==2)
  
  // **************************************************************************
  // Boundary Edges
  // **************************************************************************
  else if(tipo == -1 || tipo == 1){ // tipo = -1 or +1 (inflow or outflow)  
    iEr=0;
    iEl=0; 
    double * pdir, * sdir, pcdir[qmax];
    pdir=border.pdir;  
    if(tipo == -1){
      //printf("Inside EdgeIntegrals: tipo = -1\n");
      sdir=border.sdir;
      for(int k=0;k<qmax;k++)pcdir[k]=fluids.pressao_capilar(sdir[k]);
    }
    //else  printf("Inside EdgeIntegrals: tipo =  1\n");
    
    
    // *******************
    //     VETOR         *
    // *******************
 //   int temp;
    // *************************
    // Pressure Equation Vector
    // *************************
    
    // integral 7 

    for(int rp=0;rp<np[0];rp++){
     // e[iEr].Traco_grad_phi(a[iEr],pres,rp,gphi_r);
      e[iEr].Traco_Kgrad_phi_n(a[iEr],pres,rp,res2);
     
      prod3v(qmax,lambdat[iEr],res2,pdir,res3);
      integral(qmax,w,res3,aux);
      if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 5 " << aux << endl;
      B[e[iEr].map(pres,rp)] += aux;
    }
    
    // integral 8
    
    for(int rp=0;rp<np[0];rp++){
   //   int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);
      // if(rp_on_border){
      //e[iEr].Traco_phi_1(a[iEr],pres,pos,phi_r);
      e[iEr].Traco_phi(a[iEr],pres,rp,phi_r);
      prodvv(qmax,pdir,phi_r,res1);
      integral(qmax,w,res1,aux);
      if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 6 " << aux << endl;
      B[e[iEr].map(pres,rp)] += klt*penalty*aux; 
      
    }
    
    // inflow; acrescentar as integrais com lambda_n e p_c 
    
    if(tipo==-1){// tipo==-1
      
      // *************************
      // Pressure Equation Vector
      // *************************

      // integral 2
      
      for(int qtemp=0;qtemp<qmax;qtemp++) {
	res0[qtemp]=0.0;
	for(int ndir=0;ndir<ndim;ndir++){
	  res0[qtemp]+=Kgpc[0][ndir][qtemp]*n_e[ndir];
	}
      }
      
      //prodgn(qmax,Kgpc[0][0],Kgpc[0][1],n_e,res0);
      prodvv(qmax,lambdan[0],res0,res0);
      int temp=1;
    
      for(int rp=0;rp<np[0];rp++){
	e[iEr].Traco_phi(a[iEr],pres,rp,phi_r);
	prodvv(qmax,res0,phi_r,res1);
	integral(qmax,w,res1,aux);
	if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 7 " << aux << endl;
	B[e[iEr].map(pres,rp)] += aux*temp;
      }
      
      // integrais 4 e 9 
      
      subtvv(qmax,pc[0],pcdir,res0);// res0 eh usado em 4,6,9,10
      for(int rp=0;rp<np[0];rp++){
	//e[iEr].Traco_grad_phi(a[iEr],pres,rp,gphi_r);
	e[iEr].Traco_Kgrad_phi_n(a[iEr],pres,rp,res1);

	prod3v(qmax,lambdan[iEr],res1,res0,res3);
	integral(qmax,w,res3,aux);
	if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 8 " << aux << endl;
	B[e[iEr].map(pres,rp)] -= aux;
      }
      
      // integrais 6 e 10
      
      for(int rp=0;rp<np[0];rp++){
	e[iEr].Traco_phi(a[iEr],pres,rp,phi_r);
	prodvv(qmax,res0,phi_r,res1);
	integral(qmax,w,res1,aux);
	if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 9 " << aux << endl;
	B[e[iEr].map(pres,rp)] -= kln*penalty1*aux; // 07/08/08
      }
    }

    // OK 22/05/2008
    // ********
    // Matriz *
    // ********
    // ************************************************************************
    // Pressure equation
    // ************************************************************************
    
    for(int rp=0;rp<np[iEr];rp++){// PE_BE  rp
      e[iEr].Traco_phi(a[iEr],pres,rp,phi_r);
      //e[iEr].Traco_grad_phi(a[iEr],pres,rp,gphi_r);
      e[iEr].Traco_Kgrad_phi_n(a[iEr],pres,rp,res0);
      
      // Pressure
      for(int lp=0;lp<np[iEl];lp++){ // PE_BE_P lp
	e[iEl].Traco_phi(a[iEl],pres,lp,phi_l);
	//e[iEl].Traco_grad_phi(a[iEl],pres,lp,gphi_l);
	e[iEl].Traco_Kgrad_phi_n(a[iEl],pres,lp,res1);
	
	// integral 1

	prod3v(qmax,lambdat[iEl],res1,phi_r,res1);
	integral(qmax,w,res1,aux1);
	
	// integral 3
	
	prod3v(qmax,lambdat[iEr],res0,phi_l,res2);
	integral(qmax,w,res2,aux2);
	
	// integral 5
	
	prodvv(qmax,phi_l,phi_r,res3);
	integral(qmax,w,res3,aux3);
	
	aux = -aux1+aux2+klt*penalty*aux3;
	if(std::isnan(aux))cout << "DG_EI_IG  Float was Not a Number: aux 10 " << aux << endl;
	Ti[count]=e[iEr].map(pres,rp);
	Tj[count]=e[iEl].map(pres,lp);
	Tx[count]=aux;
	count++;
      }
      // OK 22/05/2008  
    }
  }       
  //printf(" count = %d\n",count);
};

