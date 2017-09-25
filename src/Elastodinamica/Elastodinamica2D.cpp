/*****************************************************************************/
/*****************************************************************************/
#undef __ED_run
#define __ED_run
/*****************************************************************************/
/*****************************************************************************/
# include "ED_Prob.h"
# include "PhElem.hpp"
# include "GeProb.hpp"

// Definicao da classe ED_Prob
ED_Prob::ED_Prob(Epetra_Comm& comm) : GeProb<MyElem,2,2>::GeProb(comm)
{};

ED_Prob::~ED_Prob()
{
  /*
  cout << "\nED_Prob::~ED_Prob()\n";
  // *************************************************************************
  // Liberar memoria local
  // *************************************************************************
  const int nsat=el[0].show_ptr_stdel(sat)->nn_val();// Num de modos para saturacao
  const int npres=el[0].show_ptr_stdel(pres)->nn_val();// Num de modos para pressao
  cout << "Ponto1 nsat ="<< nsat<< " npres = "<< npres << "\n";
  liberar_mem_local(nsat,npres);
  
  for(int i = 0; i < NELEM; i++) {
    el[i].finaliza_vetores();
  }
  
  
  delete [] gbtrpw;gbtrpw=nullptr;
  delete [] gbtrsn;gbtrsn=nullptr;
  
  cout << "Final de ED_Prob::~ED_Prob()\n";
  */
};


// ****************************************************************************
// ****************************************************************************
void ED_Prob::SolveElast2D_dyn(FILE * fout, double lambda, double mu)
// ****************************************************************************
{
  printf("\nSolveElast2D_dyn (COMECO)\n");
  // Compatibilidade de N_VAR (N_VAR eh uma variavel private de ED_Prob
  if(N_VAR != 2){
    printf("ERRO: N_VAR (= %d) != 2",N_VAR);
    exit(0);
  }
  // *************************************************************************
  // Parametros para a iteracao temporal
  //
  double deltaT, delta, alpha,Tf,aux0;
  double a[8];
  int prnflag=0;
  FILE *fin;
  fin=fopen("entrada","rb");
  //printf("Entre com o valor de Tf:\n");
  fscanf(fin,"%lf %*s",&Tf);
  //printf("Entre com o valor de deltaT:\n");
  fscanf(fin,"%lf %*s",&deltaT);
  //printf("Entre com o valor de delta(>= 0.5):\n");
  fscanf(fin,"%lf %*s",&delta);
  aux0 = (delta+0.5)*(delta+0.5)*0.25;
  //printf("Entre com o valor de alpha (>=%lf):\n",aux0);
  fscanf(fin,"%lf %*s",&alpha);
  fclose(fin);
  a[0]=1.0/alpha/deltaT/deltaT;
  a[1]=delta/alpha/deltaT;
  a[2]=1.0/alpha/deltaT;
  a[3]=0.5/alpha-1.0;
  a[4]=delta/alpha-1.0;
  a[5]=deltaT/2.0*(delta/alpha-2.0);
  a[6]=deltaT*(1.0-delta);
  a[7]=delta*deltaT;
  //
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  // **************************************************************************
  // INICIALIZACAO DE DADOS NOS ELEMENTOS
  // Calcula o Jacobiano e o vetor inicial de cada elemento
  // **************************************************************************
  for(int e=0;e<NELEM;e++){
    el[e].compute_J();// deve ser calculado antes do vetor

    //el[e].make_vector(0,funcao);//Calcula vetor incluindo nos internos field=0
    //el[e].make_vector(1,g);// field=1

    //printf("Termino do elemento %d\n",e);
    Stdel * pointer=el[e].show_stdel();
    int nn=pointer->nn_val();
    int NLV=el[e].show_NumLocalVars();
    int nl=nn*NLV;
    int nb=pointer->nb_val();
    int Nb=nb*NLV;
    int Ni=nl-Nb;
    // Alocacao dinamica das matrizes
    NEWMAT::Matrix * ptKcomp   = new NEWMAT::Matrix(Nb,Nb);
    NEWMAT::Matrix * ptKi_inv  = new NEWMAT::Matrix(Ni,Ni);
    NEWMAT::Matrix * ptKcKi_inv= new NEWMAT::Matrix(Nb,Ni);
    // Matrix * ptM       = new Matrix(0,nl-1,0,nl-1);
    el[e].AlocarMatrizesK(ptKcomp,ptKi_inv,ptKcKi_inv);
    // el[e].AlocarMatrizM(ptM);
  }// Fim do for sobre os elementos
  printf("Terminou Calculo dos Jacobianos. clock acumulado=%u\n",(unsigned)clock());
  // **************************************************************************
  // Termino da INICIALIZACAO dos dados dos elementos
  // **************************************************************************

  // **************************************************************************
  // Etapa inicial: Criar Matriz A
  //
  // **************************************************************************

  int nz;
  nz=el[0].show_stdel()->nb_val()*el[0].show_NumLocalVars();
  nz*=(nz*NELEM);
  printf("nz=%4d\n",nz);

  NumC=NG-NumD;// Numero de (variaveis) Conhecidas
  printf("Ponto 7: Montagem da Matriz Global\nNG =%d NumC = %d, NumD=%d\n", NG, NumC,NumD);
  // **************************************************************************
  // QUANDO TEM INCOGNITAS (NumD > 0)
  // **************************************************************************
  if(NumD>0){

    int    * Ti = new int    [nz];
    int    * Tj = new int    [nz];
    double * Tx = new double [nz];


    // ************************************************************************
    // Montagem da Matriz A em forma de Ti, Tj, Tx
    // ************************************************************************
    int NC=NumC-1;
    int ND=NumD-1;
    NEWMAT::Matrix C(NumD,NumC);
    for(int j=0;j<NumC;j++){
      for(int i=0;i<NumD;i++)C[i][j]=0.0;
    }
    int nb=ptstdel->nb_val();
    int Nb=nb*N_VAR;
    int nl=N_VAR*ptstdel->nn_val();
    int count =0;
    //printf("Ponto 7a: Loop sobre os elementos\n");
    for(int e=0;e<NELEM;e++){
      // printf("%d ",e);
      el[e].make_Kcomp_Elast2D_dyn(lambda,mu,a);
      for(int i=0;i<nb;i++){
	for(int ia=0;ia<N_VAR;ia++){
	  int ii=novoNum[el[e].map(i,ia)];
	  int il=i*N_VAR+ia;
	  if(ii>NG)
	    printf("Erro no mapeamento do elemento %d indice = %d\n",e,i);
	  if(ii<NumD){
	    for(int j=0;j<nb;j++){
	      for(int ja=0;ja<N_VAR;ja++){
		int jl=j*N_VAR+ja;
		int jj=novoNum[el[e].map(j,ja)];
		if(jj>NG)
		  printf("Erro no mapeamento do elemento %d indice = %d\n",e,j);
		if(jj<NumD){
		  Ti[count]=ii;
		  Tj[count]=jj;
		  Tx[count]=el[e].Kcomp(il,jl);
		  count++;
		}
		else  C[ii][jj-NumD]+=el[e].Kcomp(il,jl);
	      }
	    }
	  }
	}
      }
    }
    // nz=count;
    printf("\nPonto 7b: umfpack_di_triplet_to_col\nnz= %5d count = %5d\n",nz,count);
    if(nz<count) exit(0);
    nz=count;
    int      Ap[NumD+1];
    int    * Ai = new int    [nz];
    double * Ax = new double [nz];
    int    * Map= new int    [nz];
    umfpack_di_triplet_to_col(NumD,NumD,nz,Ti,Tj,Tx,Ap,Ai,Ax,Map);
    printf("ponto 7c: Terminou\nnz: inicial= %5d ; final = %5d\n",nz,Ap[NumD]);
    //for(int i=0;i<nz;i++)
    // fprintf(fout,"Ap[%d] = %10.4g\n",i,Ax[i]);

    // ************************************************************************
    // Prepara a matriz para resolucao por umfpack - Cria void * Numeric
    // ************************************************************************
    printf("Ponto 7d: umfpack_di_symbolic\n");
    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ;
    (void) umfpack_di_symbolic (NumD, NumD, Ap, Ai, Ax, &Symbolic, null, null);
    printf("Ponto 7e: umfpack_di_numeric\n");
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic);

    // ************************************************************************
    // Executar o trecho abaixo "Quando existe variavel conhecida"
    //
    //
    // ^^^^^^^^^^^Fim de NumC>0 (Trecho de variaveis conhecidas)^^^^^^^^^^^^^^^
    //} // Extra
    //} //Extra

    // 2^^^^^^^^^^^^^^^^^MONTAGEM DA MATRIZ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // ************************************************************************
    printf("Montou Matriz.\n");
    // **********************************************************************
    // ***********VETOR B****************************************************
    double B[NumD];// vetor
    Vector BSchur(0,ND);// vetor para a armazenar a contribuicao da c.d.c.
    double aux[N_VAR];
    // ************************************************************************
    // Parte iterativa comeca aqui
    // ************************************************************************
    //count=0;
    printf("\nComeco da iteracao\n");
    for(double t=0.0; t<= Tf; t+=deltaT){
      //printf("%d\n",count++);
      for(int i=0;i<NumD;i++){
 	B[i]=0.0;
 	BSchur[i]=0.0;
      }

      // condicao de contorno de forca na variavel externa 22
      //   B[novoNum[22]]=1.0;
      //   B[novoNum[2]]=1.0;
      //   B[novoNum[44]]=1.0;
      //   B[novoNum[46]]=1.0;
      //   B[novoNum[48]]=1.0;
      //   B[novoNum[50]]=1.0;

      //  for(int i=0;i<NG;i++){
      //    int ii = novoNum[i];
      //    if(ii<NumD){
      //      B[ii]=X0[i];// condicoes de contorno
      //    }
      //  }
      // Condicao de contorno inicial de solo
      if(t==0.0)B[novoNum[110*N_VAR+1]]=1.0e9;

      for(int e=0;e<NELEM;e++){
 	el[e].make_vector_Elast2D_dyn(a);
 	for(int i=0;i<nb;i++){
 	  el[e].VectorElast2D(i,aux);
 	  for(int ia=0;ia<N_VAR;ia++){
 	    int ii=novoNum[el[e].map(i,ia)];
 	    if(ii>NG)
 	      printf("Erro no mapeamento do elemento %d indice = %d\n",e,i);
 	    if(ii<NumD){
 	      B[ii]+=aux[ia];
 	    }
 	  }
 	}
      }

      if(NumC>0){
 	int NC=NumC-1;
 	// ********************************************************************
 	// Condensacao de Schur: construcao do vetor Baux
 	// Condicoes de contorno nulas X0=0
 	//printf("Montagem do vetor de Schur\n");

 	Vector Baux(0,NC);
	//double Baux[NumC];
 	for(int j=0;j<NumC;j++) Baux[j]=0.0;
 	for(int i=0;i<NG;i++){
 	  int ii=novoNum[i];
 	  if(ii>=NumD){
 	    Baux[ii-NumD]=X0[i];// deslocamentos conhecidos
 	  }
 	}
 	// Condensacao de Schur: Recalcular vetor independente B

	// for(int i=0;i<NumD;i++){
	//   aux0=0.0;
	//   for(int j=0;j<NumC;j++)aux0+=C[i][j]*Baux[j];
	//   BSchur[i] = aux0;
	// }
	BSchur = C * Baux;
 	for(int i=0;i<NumD;i++)
 	  B[i]-=BSchur[i];
 	// ^^^^^^^^^^^^^^Condensacao de Schur^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      }// <=============<

      //printf("Entrando na resolucao do sistema global\n");
      //
      // ^^^^^^^^^Fim de NumC>0 (Trecho de variaveis conhecidas)^^^^^^^^^^^^^^^
      // **********************************************************************

      // ^^^^^^^^^^^VETOR B^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      // **********************************************************************

      // **********************************************************************
      // Resolucao do sistema
      // Usando umfpack para resolver o sistema

      double x[NumD];
      double U[NG];

      (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, B, Numeric, null, null);

      //printf("Terminou resolucao do sistema global\n");
      // ^^^^^^^^^^^^^^^^SOLUCAO DO SISTEMA^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      // Resolveu sistema!
      // **********************************************************************

      // **********************************************************************
      // Passar solucao para vetor U
      for(int i=0;i<NG;i++){
 	int ii=novoNum[i];
 	if(ii<NumD){
 	  U[i]=x[ii];// Deslocamento
 	}
 	else U[i]=X0[i];
      }
      //printf("\nValores nos pontos de Gauss de X0\n");
      // ^^^^^^^^^^^^^^^^^^^ATUALIZACAO DE X0^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      // **********************************************************************
      //   printf("DESLOCAMENTO INICIAL\ncoeficientes modais\n");
      //   for(int ia=0;ia<N_VAR;ia++){
      //     printf("Campo %2d\n",ia);
      //     for(int i=0;i<NG/N_VAR;i++){
      //       int ig=i*N_VAR+ia;
      //       printf("i= %3d X0=%14.5e \n",i,X0[ig]);
      //     }
      //   }

      // **********************************************************************
      //                             OUTPUT
      // Passar solucao global para elementos calcularem valores internos
      // e fazerem output
      //fprintf(fout,"tempo = %lf\n",t);
      //printf("tempo = %lf\n",t);
      for(int e=0;e<NELEM;e++){
 	//if(e==191)prnflag=1;
 	//else prnflag=0;
 	el[e].P_eval_print_Elast2D_dyn(fout,prnflag,U,lambda,mu,a);
 	//Updates internal modes and outputs result to file
      }
      int no=112;
      fprintf(fout,"%11.4e %11.4e ",U[no*N_VAR],U[no*N_VAR+1]);
      no=114;
      fprintf(fout,"%11.4e %11.4e ",U[no*N_VAR],U[no*N_VAR+1]);
      no=116;
      fprintf(fout,"%11.4e %11.4e ",U[no*N_VAR],U[no*N_VAR+1]);
      no=118;
      fprintf(fout,"%11.4e %11.4e\n",U[no*N_VAR],U[no*N_VAR+1]);
      // ^^^^^^^^^^^^^^^^^^^^^^^^^^^OUTPUT^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      // **********************************************************************
    }
    // ^^^^^^^^^^^^^^^^TERMINOU PARTE ITERATIVA^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // ************************************************************************
    printf("\n");
    // Liberar Numeric so quando nao necessitar mais dele
    umfpack_di_free_numeric (&Numeric) ;
    // ************************************************************************

  }
  // ^^^^^^^^^^^^^^^^^^Tem incognita^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // **************************************************************************
}// FIM ///////////////////////////////////////////////////////////////////////
// ****************************************************************************

// ****************************************************************************
void PhElem::make_vector_Elast2D_dyn(const double a[])
{
  double aux;
  int NN=2*pt_stdel->nn_val();
  for(int i=0;i<NN;i++){
    aux=0.0;
    for(int j=0;j<NN;j++)
      aux+=( M[i][j] * (a[0]*u0[j]+a[2]*u1[j]+a[3]*u2[j]) );// pode usar vmath
    b0[i]=aux;
  }
  // Adiciona termo de fronteira
  int NB=2*pt_stdel->nb_val();
  for(int i=0;i<NB;i++)b0[i]+=bs[i];// pode usar vmath::vadd(NB,bs,1,b0,1,b0,1);
};


// ****************************************************************************
void PhElem::make_Kcomp_Elast2D_dyn(const double lambda, const double mu,const double a[])
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
  NEWMAT::Matrix Kb(Nb,Nb); // Mb
  NEWMAT::Matrix Ki(Ni,Ni); // Mi
  NEWMAT::Matrix Kc(Nb,Ni);
  double ** K;

 // Armazenamento dinamico de K temporario
  K = new double * [nl];
  for (int i=0; i<nl; i++){
    K[i] = new double [nl];
  }
  // Construir a matriz K
  make_K_Elast2D(lambda,mu,K);

  // Construir a matriz M
  M = new double * [nl];
  for (int i=0; i<nl; i++){
    M[i] = new double [nl];
  }
  double aux;
  for(int i=0;i<nn;i++){
    int ii=i*NumFields;
    for(int j=0;j<nn;j++){
      int jj=j*NumFields;
      aux=pt_stdel->mass(i,j,JV)*rho;
      for(int ia=0;ia<2;ia++)
	M[ii+ia][jj+ia]=aux;
      M[ii][jj+1]=0.0;
      M[ii+1][jj]=0.0;
    }
  }

  // Kb, Kc, Ki
  for(int ivar=0;ivar<2;ivar++){
    for(int n0=0;n0<nn;n0++){
      int nx=n0*2+ivar;
      for(int n1=0;n1<nn;n1++){
	for(int jvar=0;jvar<2;jvar++){
	  int ny=n1*2+jvar;
	  double aux=K[nx][ny] + a[0]*M[nx][ny]; // <===============< //
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
  for(int i=0;i<nl; i++) delete [] K[i];
  delete [] K;
}

// ****************************************************************************
// Atualizacao da velocidade e aceleracao e
// saida para arquivo dos deslocamentos
// ----------------------------------------------------------------------------
void PhElem::P_eval_print_Elast2D_dyn(FILE * fout,const int prnflag,
				    const double X[],double lambda, double mu,
				    const double a[])
{
  int nb=pt_stdel->nb_val();
  int Nb=nb*NumFields;
  int nn=pt_stdel->nn_val();
  int Nn=nn*NumFields;
  int Ni=Nn-Nb;
  int ii,gbi;
  double tu0[Nn],tu1[Nn], tu2[Nn];
  // Transfere os coeficientes dos modos de contorno para vetor local
  for(int i=0;i<nb;i++){
    for(int ia=0;ia<NumFields;ia++){
      ii=i*NumFields+ia;
      gbi=map(i,ia);
      tu0[ii]=X[gbi]*sgn[i];
    }
  }
  // Calcula os coeficientes dos modos internos
  for(int i=Nb;i<Nn;i++){
    tu0[i]=0.0;
    for(int j=0;j<Nb;j++){
      tu0[i]-=(*ptKcKi_inv)(j,i-Nb)*tu0[j];
    }
    for(int j=Nb;j<Nn;j++){
      tu0[i]+=(*ptKi_inv)(i-Nb,j-Nb)*b0[j];
    }
  }
  // Atualiza velocidade e aceleracao
  for(int i=0;i<Nn;i++){
    tu2[i]=a[0]*(tu0[i]-u0[i]) - a[2]*u1[i] - a[3]*u2[i];
    tu1[i]= u1[i] + a[6]*u2[i] + a[7]*tu2[i];
  }

  for(int i=0;i<Nn;i++){
    u0[i]=tu0[i];
    u1[i]=tu1[i];
    u2[i]=tu2[i];
  }

  if(prnflag==1){
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
  }
};
