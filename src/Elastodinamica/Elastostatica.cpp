# include "spectral.h"
// ****************************************************************************
// Elastostatica 2D
// ****************************************************************************

// ****************************************************************************
void PhElem::make_K_Elast2D(const double lambda, const double mu,double ** K)
{
  // **************************************************************************
  // OBSERVACAO: NAO MULTIPLICAR POR sgn[i]*sgn[j]; 
  // Sera feito na montagem da matriz global e na avaliacao do vetor interno
  // **************************************************************************
  int nn=pt_stdel->nn_val();
  int q0=pt_stdel->Q_val(0);
  int q1=pt_stdel->Q_val(1);
  int q2=pt_stdel->Q_val(2);
  int NGQP=q0*q1*q2;
  double Phi[NGQP][nn], phi[NGQP];
  double d0_Phi[NGQP][nn],d1_Phi[NGQP][nn];// Matrizes com as derivadas nas direcoes 0(x) e 1(y)
  double grad[NGQP][3];
  double d0_v0[NGQP], d1_v0[NGQP], d0_v1[NGQP], d1_v1[NGQP];
  // Calcular Phi e as suas derivadas d0_Phi e d1_Phi
  for(int n=0;n<nn;n++){// n=modo
    // calcular os valores de Phi_n nos pontos de quadratura
    pt_stdel->eval_Phi(n,phi);
    // Calcular as derivadas de Phi_n
    pt_stdel->Gradiente(grad,phi,ptvert,gbnmap);
    for(int m=0;m<NGQP;m++){
      Phi[m][n]=phi[m];
      d0_Phi[m][n]=grad[m][0];
      d1_Phi[m][n]=grad[m][1];
    }
  }
  
  // Construir a matriz
  double trE, T_00, T_01, T_11,wi,wj;
  int count=0;
  for(int ivar=0;ivar<2;ivar++){
    for(int n0=0;n0<nn;n0++){
      int nx=n0*2+ivar;
      for(int n1=0;n1<nn;n1++){
	for(int jvar=0;jvar<2;jvar++){
	  int ny=n1*2+jvar; 
	  for(int m=0;m<NGQP;m++){
	    if(jvar==0){
	      d0_v0[m]=d0_Phi[m][n1];
	      d1_v0[m]=d1_Phi[m][n1];
	      d0_v1[m]=0.0;
	      d1_v1[m]=0.0;
	    }
	    else{
	      d0_v0[m]=0.0;
	      d1_v0[m]=0.0;
	      d0_v1[m]=d0_Phi[m][n1];
	      d1_v1[m]=d1_Phi[m][n1];	
	    }
	  }
	  // Integrar sobre o elemento
	  int m=0;// indice do ponto de integracao de Gauss
	  double aux=0.0;
	  for(int j=0;j<q1;j++){
	    wj=pt_stdel->w_val(1,j);
	    double aaux=0.0;
	    for(int i=0;i<q0;i++){// i runs fastest
	      wi=pt_stdel->w_val(0,i);
	      // traco de E
	      trE=d0_v0[m]+d1_v1[m];
	      // Tensoes
	      T_00=lambda*trE+2.0*mu*d0_v0[m];
	      T_11=lambda*trE+2.0*mu*d1_v1[m];
	      T_01=mu*(d0_v1[m]+d1_v0[m]);
	      if(ivar==0)
		aaux+=(T_00*d0_Phi[m][n0]+T_01*d1_Phi[m][n0])*wi*JV[m];//<---JV---<
	      else
		aaux+=(T_01*d0_Phi[m][n0]+T_11*d1_Phi[m][n0])*wi*JV[m];//<---JV---<
	      m++;// acrescenta o contador do ponto gaussiano
	    }
	    aux+=(wj*aaux);
	  } 
	  K[nx][ny]=aux;
	}
      }
    }
  }
};

// ***************************************************************************
void PhElem::P_eval_print_Elast2D(FILE * fout, 
				const double X[],double lambda, double mu)
{
  int nb=pt_stdel->nb_val();
  int Nb=nb*NumFields;
  int nn=pt_stdel->nn_val();
  int Nn=nn*NumFields;
  int Ni=Nn-Nb;
  int ii,gbi;
  // Transfere os coeficientes dos modos de contorno para vetor local
  for(int i=0;i<nb;i++){
    for(int ia=0;ia<NumFields;ia++){
      ii=i*NumFields+ia;
      gbi=map(i,ia);
      u0[ii]=X[gbi]*sgn[i];
    }
  }
  // Calcula os coeficientes dos modos internos
  for(int i=Nb;i<Nn;i++){
    u0[i]=0.0;
    for(int j=0;j<Nb;j++){
      u0[i]-=(*ptKcKi_inv)(j,i-Nb)*u0[j];
    }
    for(int j=Nb;j<Nn;j++){
      u0[i]+=(*ptKi_inv)(i-Nb,j-Nb)*b0[j];
    }
  }
  
  double utemp[nn];
  int D=pt_stdel->Q_val(0)*pt_stdel->Q_val(1)*pt_stdel->Q_val(2);
  double x[D],y[D],z[D],ftemp[D];
  double f[D][NumFields]; 
  
  for(int numf=0;numf<NumFields;numf++){
    for(int i=0;i<nn;i++){
      ii=i*NumFields+numf;
      utemp[i]=u0[ii];
    }
    pt_stdel->printtoarray(utemp,ptvert,gbnmap,x,y,z,ftemp);
    for(int i=0;i<D;i++){
      f[i][numf]=ftemp[i];
    }
  }
  
  for(int i=0;i<D;i++)
    fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x[i],y[i],z[i],f[i][0],f[i][1]); 
  
};
// ****************************************************************************
// ****************************************************************************
void PhElem::VectorElast2D(int i,double sum[])
{
  // *********************************************************
  // i = numero local do no; cada no tem NumFields variaveis *
  // *********************************************************
  
  int nb=pt_stdel->nb_val();// numero de nos de contorno
  int ni=pt_stdel->nn_val()-nb;// numero de nos internos
  int Ni=ni*NumFields;
  int Nb=nb*NumFields;  
  int ia;
  for(ia=0;ia<NumFields;ia++){
    int ii = NumFields*i + ia;
    sum[ia]=b0[ii];
    for(int j=0; j<Ni; j++){
      sum[ia]-=(*ptKcKi_inv)(ii,j)*b0[j+Nb];
    }
  }
  for(ia=0;ia<NumFields;ia++)
    sum[ia]*=sgn[i];
};
// ****************************************************************************
// ****************************************************************************
void PhElem::VectorElast2D(double vector[])
{
  // *********************************************************
  // i = numero local do no; cada no tem NumFields variaveis *
  // *********************************************************

  int nb=pt_stdel->nb_val();// numero de nos de contorno
  int ni=pt_stdel->nn_val()-nb;// numero de nos internos
  int Ni=ni*NumFields;
  int Nb=nb*NumFields;
  int ia;
  for(int i=0;i<Nb;i++){
    vector[i]=b0[i];
    for(int j=0; j<Ni; j++){
      vector[i]-=(*ptKcKi_inv)(i,j)*b0[j+Nb];
    }
  }
  for(int i=0;i<nb;i++)
    for(ia=0;ia<NumFields;ia++)
      vector[i*NumFields+ia]*=sgn[i];
};
// ****************************************************************************
double PhElem::Kcomp(const int i, const int j)
{
  return ((*ptKcomp)(i,j)*sgn[(i/NumFields)]*sgn[(j/NumFields)]);
};

// ****************************************************************************
void PhElem::teste()
{
  NumFields=0;
  pt_stdel->teste(NumFields);
  printf("NumFields=%d\n",NumFields);
};

// ****************************************************************************
void PhElem::make_Kcomp_Elast2D(const double lambda, const double mu)
{
  // **************************************************************************
  // OBSERVACAO: NAO MULTIPLICAR POR sgn[i]*sgn[j]; 
  // Sera feito na montagem da matriz global e na avaliacao do vetor interno
  // **************************************************************************
  int Nb=2*pt_stdel->nb_val();
  int nn=pt_stdel->nn_val();
  int Ni=2*nn-Nb;
  int nl = 2*nn;
  // Matrizes locais
  Matrix Kb(0,Nb-1,0,Nb-1); // Mb
  Matrix Ki(0,Ni-1,0,Ni-1); // Mi
  Matrix Kc(0,Nb-1,0,Ni-1);
  double ** K;

 // Armazenamento dinamico de K temporario
  K = new double * [nl];
  for (int i=0; i<nl; i++){
    K[i] = new double [nl];
  }
  // Construir a matriz K
  make_K_Elast2D(lambda,mu,K);

  // Kb, Kc, Ki
  for(int ivar=0;ivar<2;ivar++){
    for(int n0=0;n0<nn;n0++){
      int nx=n0*2+ivar;
      for(int n1=0;n1<nn;n1++){
	for(int jvar=0;jvar<2;jvar++){
	  int ny=n1*2+jvar; 
	  double aux=K[nx][ny] ; // <===============< //
	  if(nx<Nb){
	    if(ny<Nb) Kb(nx,ny)=aux;
	    else Kc(nx,ny-Nb)=aux;
	  }
	  else if(ny>=Nb) Ki(nx-Nb,ny-Nb)=aux;
	}
      }
    }
  }
  *ptKi_inv=Ki.Inverse();
  *ptKcKi_inv=Kc*(*ptKi_inv);
  *ptKcomp=Kb - ( ( Kc * (*ptKi_inv) ) * Kc.Transpose());

  // Free memory of K
  for(int i=0;i<nl; i++) delete K[i];
  delete K;

};

// ****************************************************************************
void GeProb::SolveElast2D(FILE * fout, double lambda, double mu)
// ****************************************************************************
{  
  // **************************************************************************
  // Assembly of vector and matrices
  // **************************************************************************
  // **************************************************************************
  // INICIALIZACAO DE DADOS NOS ELEMENTOS
  // Calcula o Jacobiano e o vetor inicial de cada elemento 
  // **************************************************************************
  for(int e=0;e<NELEM;e++){
    el[e].compute_J();// deve ser calculado antes do vetor
    //printf("calculou J (Jacobiano)\n");
        
    //el[e].make_vector(0,funcao);//Calcula vetor incluindo nos internos field=0
    //el[e].make_vector(1,g);// field=1
    
    //printf("Termino do elemento %d\n",e);
    int nn=el[e].show_stdel()->nn_val();
    int nl=2*nn;
    int nb=el[e].show_stdel()->nb_val();
    int Nb=2*nb;
    int Ni=nl-Nb;
    // Alocacao dinamica das matrizes
    Matrix * ptKcomp   = new Matrix(0,Nb-1,0,Nb-1);
    Matrix * ptKi_inv  = new Matrix(0,Ni-1,0,Ni-1);
    Matrix * ptKcKi_inv= new Matrix(0,Nb-1,0,Ni-1);
    el[e].AlocarMatrizesK(ptKcomp,ptKi_inv,ptKcKi_inv);
  }// Fim do for sobre os elementos
  printf("Terminou Calculo dos Jacobianos. clock acumulado=%u\n",(unsigned)clock()); 
  // **************************************************************************
  // Etapa inicial: Criar Matriz A e vetores B0 e B1 com condicoes de contorno
  //                e condicoes iniciais nulas.
  // **************************************************************************
  
  int nz;
  nz=el[0].show_stdel()->nb_val()*el[0].show_NumFields();
  nz*=(nz*NELEM);
  printf("nz=%4d\n",nz);

  NumC=NG-NumD;// Numero de (variaveis) Conhecidas
  printf("ponto 7 NG =%d NumC = %d, NumD=%d\n", NG, NumC,NumD);
  if(NumD>0){// Tem incognita
    int ND=NumD-1;
    double B[NumD];// vetor
    Vector BSchur(0,ND);// vetor para a armazenar a contribuicao da c.d.c.
    int NC=NumC-1;
    Matrix C(0,ND,0,NC);
    Vector Baux(0,NC);
    for(int j=0;j<=NC;j++){
      Baux[j]=0.0;
      for(int i=0;i<=ND;i++)C(i,j)=0.0;
    } 
    int    * Ti = new int    [nz];
    int    * Tj = new int    [nz];
    double * Tx = new double [nz];
   
    // ***************************1********************************************
    for(int i=0;i<NumD;i++){
      B[i]=0.0;
      BSchur[i]=0.0;
    }
    
    // ***************************1********************************************
    // condicao de contorno de forca na variavel externa 22 
		//   B[novoNum[22]]=1.0;
		//   B[novoNum[2]]=1.0;
		//   B[novoNum[44]]=1.0;
   //   B[novoNum[46]]=1.0; 
   //   B[novoNum[48]]=1.0; 
   //   B[novoNum[50]]=1.0;
   
   //    for(int i=0;i<NG;i++){
   //      int ii = novoNum[i];
   //      if(ii<NumD){
   //        B[ii]=X0[i];
   //        //Bsurface[ii]=X0[i];
   //        fprintf(fout3,"Valor inicial do B[%2d] = %11.4e\n",ii,X0[i]);
   //      }
   //    }
    B[novoNum[110*NFIELDS+1]]=1.0e2;//Condicao de contorno para a malha solo.neu
    // 2***********************************************************************
    int nb=ptstdel->nb_val();
    int Nb=nb*NFIELDS;
    double aux[NFIELDS];
    int nl=NFIELDS*ptstdel->nn_val();
    int count =0;
    for(int e=0;e<NELEM;e++){
      //printf("Montagem da matriz do Elemento %d\n",e);
      el[e].make_Kcomp_Elast2D(lambda,mu);
      // Montagem do vetor b0 ( = b0+bs)
      for(int ia=0;ia<NFIELDS;ia++)
	el[e].make_vector(ia,force0);
      for(int i=0;i<nb;i++){ 
	el[e].VectorElast2D(i,aux);
	for(int ia=0;ia<NFIELDS;ia++){
	  int ii=novoNum[el[e].map(i,ia)];
	  int il=i*NFIELDS+ia;
	  if(ii>NG)
	    printf("Erro no mapeamento do elemento %d indice = %d\n",e,i);
	  if(ii<NumD){
	    B[ii]+=aux[ia];
	    for(int j=0;j<nb;j++){
	      for(int ja=0;ja<NFIELDS;ja++){
		int jl=j*NFIELDS+ja;
		int jj=novoNum[el[e].map(j,ja)];
		if(jj>NG) 
		  printf("Erro no mapeamento do elemento %d indice = %d\n",e,j);
		if(jj<NumD){
		  Ti[count]=ii;
		  Tj[count]=jj;
		  Tx[count]=el[e].Kcomp(il,jl);
		  count++;
		}
		else  C(ii,jj-NumD)+=el[e].Kcomp(il,jl);
	      }
	    }
	  }
	}
      }
    }
    // nz=count;
    printf("ponto 7b nz= %5d count = %5d\n",nz,count);
    for(int s=0;s<NumD;s++)if(B[s]!=0.0)fprintf(fout3,"B[%d]=%17.10e\n",s,B[s]);
    if(nz<count) exit(0);
    nz=count;
    int      Ap[NumD+1];
    int    * Ai = new int    [nz];
    double * Ax = new double [nz];
    int    * Map= new int    [nz];
    umfpack_di_triplet_to_col(NumD,NumD,nz,Ti,Tj,Tx,Ap,Ai,Ax,Map);
    printf("ponto 7c nz: inicial= %5d ; final = %5d\n",nz,Ap[NumD]);
    // for(int i=0;i<nz;i++)
    //  fprintf(fout,"Ax[%d] = %10.4g\n",i,Ax[i]);
    
    if(NumC>0){
      // **********************************************************************
      // construcao do vetor Baux
      // Condicoes de contorno nulas X0=X1=0
      printf("Montagem do vetor de Schur\n");
      for(int i=0;i<NG;i++){
	int ii=novoNum[i];
	if(ii>=NumD){
	  Baux(ii-NumD)=X0[i];
	}
      }
      // Condensacao de Schur
      BSchur = C * Baux;
      for(int i=0;i<NumD;i++)
	B[i]-=BSchur[i];
      // **********************************************************************
    }// Fim de NumC>0
    printf("Terminou calculo da matriz global\nEntrando na resolucao do sistema global\n");
    // Usando umfpack para resolver o sistema
    double x[NumD];
    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ; printf("ponto 7d\n");
    (void) umfpack_di_symbolic (NumD, NumD, Ap, Ai, Ax, &Symbolic, null, null);
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic);
    (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, B, Numeric, null, null);
    umfpack_di_free_numeric (&Numeric) ;
    printf("Terminou resolucao do sistema global\n");
    for(int i=0;i<NG;i++){
      int ii=novoNum[i];
      if(ii<NumD){
	X0[i]=x[ii];// Deslocamento
      }
    }
   //   printf("DESLOCAMENTO INICIAL\ncoeficientes modais\n");
   //   for(int ia=0;ia<NFIELDS;ia++){
   //     printf("Campo %2d\n",ia);
   //     for(int i=0;i<NG/NFIELDS;i++){
   //       int ig=i*NFIELDS+ia;
   //       printf("i= %3d X0=%14.5e \n",i,X0[ig]);
   //     }
   //   }
    printf("\nValores nos pontos de Gauss de X0\n");
    for(int e=0;e<NELEM;e++){
      el[e].P_eval_print_Elast2D(fout,X0,lambda,mu);//Outputs result to file
    }
    fprintf(fout,"\n\n");
    
  }
  // ^^^^^^^^^^^^^^^^^^Tem incognita^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
}
// ****************************************************************************
