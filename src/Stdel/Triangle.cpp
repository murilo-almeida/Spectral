#include "spectral.h"
//#define PRINTF_ON
// ****************************************************************************
// Class Triangle
// ****************************************************************************
// ****************************************************************************
Triangle::Triangle(int p0,int q0)
{
  //  printf("Iniciando Triangle\n");
  tipo=2;
  ndim=2;
  nv=3;
  ne=3;
  nf=0;
  nborder=3;
  vtk_type=5;
  emapi = new int [ne+1];
  emapv = new int [ne*(p0+1)];
  bmapi = emapi;
  bmapv = emapv;
  set(p0, q0);
};
// ****************************************************************************
Triangle::~Triangle()
{
  //  printf("Destruindo Triangle\n");
  delete [] emapi; emapi=nullptr;
  delete [] emapv; emapv=nullptr;
  bmapi=nullptr;
  bmapv=nullptr;
  delete [] D_Phi_val; D_Phi_val=nullptr;
  
  //libera memoria de ind_mode_
  for (int i=0;i<=P[0]+1;++i){
    delete [] ind_mode_[i];  ind_mode_[i] = nullptr;
  }
  delete [] ind_mode_;  ind_mode_ = nullptr;
  // Stdel::~Stdel();
};
// ****************************************************************************

// Inicializa os dados do elemento
void Triangle::set(int p0, int q0)
{
  int i,j,k;
  int a=3;
  // Mapeamento das arestas
  // Indices que marcam o inicio das arestas
  emapi[0]=0;
  emapi[1]=p0+1;
  emapi[2]=2*(p0+1);
  emapi[3]=3*(p0+1);

  //printf("Em Triangle::set(%d,%d)\n",p0,q0);
  P[0]=p0; Q[0]=q0; gqt[0]=3;//Gauss_Lobatto_Jacobi (inclui os dois vertices)
  P[1]=p0; Q[1]=q0; gqt[1]=2;//Gauss_Radau_Jacobi (inclui so o primeiro vertice)
  P[2]=0;  Q[2]=1;  gqt[2]=1;// Terceira dimensao
  NGQP=Q[0]*Q[1]*Q[2];
  qborder = q0;
  
  //aloca memoria para ind_mode_
  ind_mode_ = new int * [P[0]+2];
  for (int i=0;i<=P[0]+1;++i){
    ind_mode_[i] = new int [P[1]+1];
  }
  
#ifdef PRINTF_ON
  printf("Em Triangle::set apos Sair de Triangle::setStdel\n");
  printf("P: %d %d\nQ: %d %d\ngqt : %d %d\n",P[0],P[1],Q[0],Q[1],gqt[0],gqt[1]);
#endif

  // Mapeamento de nos, modos e arestas
  mode_[0].set_mode(0,0,0); ind_mode_[0][0]=0;// A
  mode_[1].set_mode(P[0],0,0); ind_mode_[P[0]][0]=1;// B
  mode_[2].set_mode(P[0]+1,P[1],0); ind_mode_[P[0]+1][P[1]]=2;// C  Ponto Colapsado

  // AB
  k=0;
  emapv[k]=0;// Vertice A
  k++;
  for(i=1;i<P[0];i++){
    mode_[a].set_mode(i,0,0); ind_mode_[i][0]=a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=1;// Vertice B

  // BC
  k=emapi[1];
  emapv[k]=1; // Vertice B
  k++;
  for(j=1;j<P[1];j++){
    mode_[a].set_mode(P[0],j,0);ind_mode_[P[0]][j]=a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=2;// Vertice C

  // AC
  k=emapi[2];
  emapv[k]=0; // Vertice A
  k++;
  for(j=1;j<P[1];j++){
    mode_[a].set_mode(0,j,0); ind_mode_[0][j]=a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=2;// Vertice C

  nb=a;
  // Interior modes : j runs fastest
  for(i=1; i<P[0]; i++){
    for(j=1; j<P[1]-i;j++){
      mode_[a].set_mode(i,j,0);ind_mode_[i][j]=a;
      a++;
    }
  }
  nn=a;
  // if(nn==nmode)printf("OK! Calculo do nmode correto: %d, nn=%d\n",nmode,nn);
  if(nb>MAXNB ||nn>MAXMODES  ||nn-nb>MAXNI){
    printf("Ajuste parametros no spectral.h\n");
    printf("Numero de nos de fronteira do elemento (nb = %d)  (MAXNB = %d)\n",nb,MAXNB);
    printf("Numero de nos do elemento nn(= %d)                (MAXMODES  = %d)\n",nn,MAXMODES);
    printf("Numero de nos internos do elemento (ni = %d)      (MAXNI = %d)\n",nn-nb,MAXNI);
    exit(0);
  }
  // Calcula os parametros de Gauss
  gauss_parameters_default(); //alpha=beta=0.0
#ifdef PRINTF_ON
  printf("passou gauss_parameters\n");
  printf("Saindo de Triangle::set: nv = %d  nn = %d (< MAXMODES = %d) nb = %d (< MAXNB = %d)\n",nv,nn,MAXMODES,nb,MAXNB);
	printf("Triangle::set nn = %d nb= %d\n",nn,nb);
  
#endif
 
  // Construcao da matriz Phi_val[nn][Q[0]*Q[1]]
  double eta1,eta2;
  double Fa[Q[0]],Fb[Q[1]];
  for(int m=0;m<nn;m++){
   
    int p= mode_[m].p_val();
    int q= mode_[m].q_val();
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      Fb[j]=Psib(P[0],P[1],p,q,eta2);
    } 
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      Fa[i]=Psia(P[0],p,eta1);
    }

    int n=0;
    for(j=0;j<Q[1];j++){
      for(i=0;i<Q[0];i++){// i runs fastest
				Phi_val[m][n++]=Fa[i]*Fb[j];
      }
    }
  }
};

// ****************************************************************************
void Triangle::print(FILE * fout)
{
  Stdel::printStdel();
  fprintf(fout,"Triangle::print\n");
  for(int i=0;i<nn_val(); i++){
    fprintf(fout,"no= %4d ",i);
    mode_[i].print(fout);
  }
  fprintf(fout,"\n");
};

// *******************************************************************
// Mass Matrix entry M(m1, m2)
// *******************************************************************
// Evaluates the inner product of Phi_m1 by Phi_m2

double Triangle::mass(int m1, int m2, const double JV[])
{
  int p,q,r,s;
 // int flag1, flag2;
  int j, i;
  double aux, aaux, Fa,Fb,Ga,Gb;
  // printf("\nEntrou Triangle::mass(%d, %d)\n",m1,m2);
  // indexes of point m1
  p= mode_[m1].p_val();
  q= mode_[m1].q_val();
  // indexes of point m2
  r= mode_[m2].p_val();
  s= mode_[m2].q_val();
  aux=0.0;
  for(j=0;j<Q[1];j++){
    aaux=0.0;
    for(i=0;i<Q[0];i++){
      Fa=Psia(P[0], p,xGQ[0][i]);//primeiro
      Ga=Psia(P[0], r,xGQ[0][i]);
      aaux+=wGQ[0][i]*Fa*Ga;//*JV[i+Q[0]*j];// <==multiplica pelo Jacobiano-< 
    }
    Fb=Psib(P[0], P[1], p, q, xGQ[1][j]);
    Gb=Psib(P[0], P[1], r, s, xGQ[1][j]);
    aux+=wGQ[1][j]*Fb*Gb*aaux;
  }
  //printf("Saiu Triangle::mass(%d, %d,JV)\n",m1,m2);
  return (aux*JV[0]);// <----------multiplica pelo Jacobiano----<
}
// ****************************************************************************
// Creates the local matrices using static condensation
// It uses NEWMAT matrices and operators
// ****************************************************************************
void Triangle::make_local_matrices()
{
  int ni;
  int i,ii,j,jj;
  ni=nn-nb;
#ifdef PRINTF_ON
  printf("Making local matrices 1: nb = %d  ni= %d\n", nb,ni);
#endif

  double M[nn][nn];
  int p,q,r,s;
 // int flag1, flag2;
  double aux, aaux, Fa,Fb,Ga,Gb;
  int m1,m2;
  for(m1=0;m1<nn;m1++){
    // indexes of point m1
    p= mode_[m1].p_val();
    q= mode_[m1].q_val();
    for(m2=m1;m2<nn;m2++){
      // indexes of point m2
      r= mode_[m2].p_val();
      s= mode_[m2].q_val();
      aux=0.0;
      for(j=0;j<Q[1];j++){
        aaux=0.0;
        for(i=0;i<Q[0];i++){
          Fa=Psia(P[0], p,xGQ[0][i]);
          Ga=Psia(P[0], r,xGQ[0][i]);
          aaux+=wGQ[0][i]*Fa*Ga;
        }
        Fb=Psib(P[0], P[1], p, q, xGQ[1][j]);
        Gb=Psib(P[0], P[1], r, s, xGQ[1][j]);
        aux+=wGQ[1][j]*Fb*Gb*aaux;
      }
      M[m1][m2]=aux;
      M[m2][m1]=aux;
    }
  }
  if(ni==0) {
    for(int i=0;i<nb;i++){
      for(int j=0;j<nb;j++)
        Mb_comp[i][j]=M[i][j];
    }
  }
  else {
  #ifdef _NEWMAT
    // newmat
    NEWMAT::Matrix mb(nb,nb);
    NEWMAT::Matrix mb_c(nb,nb);
    NEWMAT::Matrix mi(ni,ni);
    NEWMAT::Matrix mc(nb,ni);
    NEWMAT::Matrix mi_inv(ni,ni);
    NEWMAT::Matrix mcmi_inv(nb,ni);
  #endif

    // Matriz Mb
    for(i=0;i<nb;i++){
      for(j=0;j<nb;j++) mb.element(i,j)=M[i][j];
      // Matrix Mc
      for(j=nb;j<nn;j++){
        jj=j-nb;
        mc.element(i,jj)=M[i][j];
      }
    }
    // Matrix Mi
    for(i=nb;i<nn;i++){
      ii=i-nb;
      mi.element(ii,ii)=M[i][i];
      for(j=i+1;j<nn;j++){
        jj=j-nb;
        mi.element(ii,jj)=M[i][j];
        mi.element(jj,ii)=M[i][j];
      }
    }

  #ifdef _NEWMAT
    mi_inv=mi.i(); // newmat
    mcmi_inv = mc * mi_inv;
    mb_c = mb - (mcmi_inv * mc.t()); // newmat
  #endif

    for(int i=0;i<nb;i++){
      for(int j=0;j<nb;j++)
        Mb_comp[i][j]=mb_c.element(i,j);
      for(int j=nb;j<nn;j++){
        jj=j-nb;
        McMi_inv[i][jj]=mcmi_inv.element(i,jj);
      }
    }
    //printf("Making local matrices 3\n");
    for(int i=nb;i<nn;i++){
      ii=i-nb;
      for(int j=nb;j<nn;j++){
        jj=j-nb;
        Mi_inv[ii][jj]=mi_inv.element(ii,jj);
      }
    }
  }
  #ifdef PRINTF_ON
    printf("Saindo Triangle::make_local_matrices: nb = %d  ni= %d\n", nb,ni);
  #endif
};

// ****************************************************************************

void Triangle::make_mass_matrices(int NFields)
{/*
  //printf("Triangle make_mass_matrices\n");
  // ***************************
  make_local_matrices();
  // ***************************
  if(NFields > 1)duplicar_mass_matrix(NFields);
#ifdef PRINTF_ON
  printf("Saindo Triangle::make_mass_matrices");
  printf("(%2d): nb = %d  ni= %d\n",NFields,nb,nn-nb);
#endif
  */
};

// ****************************************************************************
// Sets Gauss parameters for the triangle
// Ja em coordenadas eta1 e eta2, incluindo o termo (1-eta2)/2
// (Gauss-Jacobi quadratura) Ver Apendice B de Karniadakis & Sherwin
// ****************************************************************************
void Triangle::gauss_parameters_default()
{
  double x[MAXQ];
  double wtemp[MAXQ];
  double Dtemp[MAXQ][MAXQ];

  // Valores para a terceira dimensao: 1 ponto (eta3=-1.0) com w=1.0
  xGQ[2][0]=-1.0;
  wGQ[2][0]=1.0;
  // Valores para as duas dimensoes
  for(int i=0;i<ndim;i++){
    if(gqt[i]==1)
      {
				//gqt[i]=1 Gauss-Jacobi, so pontos interiores.
				Gauss_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
      }
			else if(gqt[i]==2){
				// Ja em coordenadas eta1 e eta2, incluindo o termo (1-eta2)/2
				//             alpha=1.0, beta=0.0\/ \/\/\/\/****************************
				Gauss_Radau_Jacobi_parameters(Q[i],1.0,0.0,x,wtemp,Dtemp);
				// ********************************^^^^^^^^^*****************************
				for(int j=0; j<Q[i]; j++){
					xGQ[i][j]=x[j];
					// *****\/\/\/\/\/\/***************************************************
					wGQ[i][j]=wtemp[j]/2.0;
					// *****^^^^^^^^^^^^***************************************************
				}
			}
			else{//gqt[i]=3
				Gauss_Lobatto_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
				for(int j=0; j<Q[i]; j++){
					xGQ[i][j]=x[j];
					wGQ[i][j]=wtemp[j];
      }
    }
    for(int k=0;k<Q[i];k++)
      for(int l=0;l<Q[i];l++){
				D[k][l][i]=Dtemp[k][l];
      }
		}
};

// ****************************************************************************
// OBSERVACAO: MULTIPLICAR POR JV
// ****************************************************************************
void Triangle::vector_of_integral_of_f_Phi_dv(double vec[],
				    double (*func)(double, double, double), 
				    const Vertice vert[], const int map[],
				    const double JV[])
{
  int n,p,q;
  int j, i;
  double aux, Fa,Fb,eta1,eta2,x1,x2,x3;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  int Pdim=P[0]+2;// inclui ponto colapsado (caso especial)
  double fp_xi2[Pdim][Q[1]];
  
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
   
  // construir matriz temporaria fp_xi2
  for(p=0;p<Pdim;p++){
    // fazer loop sobre xi2
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      eaux=(1.0-eta2)/4.0;
      haux=(1.0+eta2)/2.0;
      //fazer soma em xi1
      aux=0.0;
      for(i=0;i<Q[0];i++){
        // coordenadas
        eta1=xGQ[0][i];
        x1=eaux*(faux1-eta1*gaux1)+haux*xc;
        x2=eaux*(faux2-eta1*gaux2)+haux*yc;
        x3=eaux*(faux3-eta1*gaux3)+haux*zc;
        Fa=Psia(P[0], p,eta1);
        aux+=(wGQ[0][i]*Fa*func(x1,x2,x3))*JV[i+Q[0]*j];// <======Jacobiano==<
      }
      fp_xi2[p][j]=aux;
    }
  }
  
  // Multiplicar o vetor pela matriz temporaria
  for(n=0;n<nn;n++){
    p= mode_[n].p_val();
    q= mode_[n].q_val();
    aux=0.0;
    // fazer soma sobre xi2
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      Fb=Psib(P[0], P[1], p, q, eta2);
      aux+=wGQ[1][j]*Fb*fp_xi2[p][j];
    }
    vec[n]=aux;
  }
  
  // ------Multiplicar pelo Jacobiano----
  //for(n=0;n<nn;n++)vec[n]*=JV[0];// <------Multiplicar pelo Jacobiano----<
};
// ****************************************************************************
void Triangle::vector_of_integral_of_f_Phi_dv(double vec[],
				    const double func[],
				    const double JV[])
{
  int n,p,q;
  int j, i;
  double aux, Fa,Fb,eta1,eta2;//x1,x2,x3;
//  double eaux,haux;//faux1,faux2,faux3,gaux1,gaux2,gaux3,
 // double xa,ya,za,xb,yb,zb;//xc,yc,zc;
  int Pdim=P[0]+2;// inclui ponto colapsado (caso especial)
  double fp_xi2[Pdim][Q[1]];
  
   
  // construir matriz temporaria fp_xi2
  for(p=0;p<Pdim;p++){
    // fazer loop sobre xi2
    int count=0;
    for(j=0;j<Q[1];j++){
  
      //fazer soma em xi1
      aux=0.0;
      for(i=0;i<Q[0];i++){
	// coordenadas
        eta1=xGQ[0][i];
    
        Fa=Psia(P[0], p,eta1);
        aux+=(wGQ[0][i]*Fa*func[count])*JV[i+Q[0]*j];// <===Jacobiano==<
        count++;
      }
      fp_xi2[p][j]=aux;
    }
  }
  
  for(n=0;n<nn;n++){
    p= mode_[n].p_val();
    q= mode_[n].q_val();
    aux=0.0;
    // fazer soma sobre xi2
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      Fb=Psib(P[0], P[1], p, q, eta2);
      aux+=wGQ[1][j]*Fb*fp_xi2[p][j];
    }
    vec[n]=aux;
  }
  
  // ------Multiplicar pelo Jacobiano----
  //for(n=0;n<nn;n++)vec[n]*=JV[0];// <------Multiplicar pelo Jacobiano----<
};
// ****************************************************************************

// ****************************************************************************
void Triangle::printtofile(FILE * fout,const double u[],
			   double (*func)(double,double,double), 
			   const Vertice vert[], const int map[])
{
  double aux,eta1,eta2,x1,x2,x3;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  double ftemp[Q[0]*Q[1]];
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
  int m1=0;
  evalGQ(ftemp,u);
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/4.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
      eta1=xGQ[0][i];
      x1=eaux*(faux1-eta1*gaux1)+haux*xc;
      x2=eaux*(faux2-eta1*gaux2)+haux*yc;
      x3=eaux*(faux3-eta1*gaux3)+haux*zc;
      aux=ftemp[m1++];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
    }
  }
};
// ****************************************************************************
void Triangle::printtofile(FILE * fout,const double u[], 
			   const Vertice vert[], const int map[])
{
  double aux,eta1,eta2,x1,x2,x3;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  double ftemp[Q[0]*Q[1]];
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
  int m1=0;
  evalGQ(ftemp,u);
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/4.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
      eta1=xGQ[0][i];
      x1=eaux*(faux1-eta1*gaux1)+haux*xc;
      x2=eaux*(faux2-eta1*gaux2)+haux*yc;
      x3=eaux*(faux3-eta1*gaux3)+haux*zc;
      aux=ftemp[m1++];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
    }
  }
};

// ****************************************************************************
void Triangle::printGQtofile(FILE * fout,const double ftemp[],
			     const double ftemp1[],
			     const Vertice vert[], const int map[])
{
  double aux, aux1,eta1,eta2,x1,x2,x3;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
  int m1=0;
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/4.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
      eta1=xGQ[0][i];
      x1=eaux*(faux1-eta1*gaux1)+haux*xc;
      x2=eaux*(faux2-eta1*gaux2)+haux*yc;
      x3=eaux*(faux3-eta1*gaux3)+haux*zc;
      aux=ftemp[m1];
      aux1=ftemp1[m1++];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,aux1);
    }
  }
};
// ****************************************************************************
// ****************************************************************************
void Triangle::printwGQtofile(FILE * fout,
			      const Vertice vert[],
			      const int map[],
			      const double JV[])
{
  double aux,eta1,eta2,x1,x2,x3;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
  int m1=0;
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/4.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
      eta1=xGQ[0][i];
      x1=eaux*(faux1-eta1*gaux1)+haux*xc;
      x2=eaux*(faux2-eta1*gaux2)+haux*yc;
      x3=eaux*(faux3-eta1*gaux3)+haux*zc;
      aux=wGQ[0][i]*wGQ[1][j]*JV[m1++];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
    }
  }
};
// ****************************************************************************
void Triangle::printtoarray(const double u[], 
			    const Vertice vert[], const int map[],
			    double x[], double y[], double z[], double ftemp[])
{
  double eta1,eta2;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
  int m1=0;
  evalGQ(ftemp,u);
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/4.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
      eta1=xGQ[0][i];
      /*x1*/x[m1]=eaux*(faux1-eta1*gaux1)+haux*xc;
      /*x2*/y[m1]=eaux*(faux2-eta1*gaux2)+haux*yc;
      /*x3*/z[m1]=eaux*(faux3-eta1*gaux3)+haux*zc;
      //aux=ftemp[m1];
      m1++;
      //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
    }
  }
}
// ****************************************************************************
// Calcular os valores da funcao nos pontos de Gauss
// ****************************************************************************
void Triangle::evalGQ(double f0[],double f1[],
			const double uh0[],const double uh1[])
{
  double aux0,aux1, Fa,Fb,eta1,eta2;
  int a;
  int i,j,p,q,k;
  int n=0;
  int Pdim=P[0]+2;	
  double ftemp0[Pdim],ftemp1[Pdim];
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    // ************************************************************************
    // Construcao dos vetores temporarios
    // ************************************************************************
    for(k=0;k<Pdim;k++){
      ftemp0[k]=0.0;
      ftemp1[k]=0.0;
    }
    // A
    a=0;p=0;q=0;
    Fb=Psib(P[0],P[1],p,q,eta2);
    ftemp0[p]+=Fb*uh0[a];
    ftemp1[p]+=Fb*uh1[a];
    // B
    a=1;p=P[0];q=0;
    Fb=Psib(P[0],P[1],p,q,eta2);
    ftemp0[p]+=Fb*uh0[a];
    ftemp1[p]+=Fb*uh1[a];
    // C Ponto colapsado: p=P[0]+1 !!!!!
    a=2;
    p=P[0]+1;q=P[1]; // Ponto com p=P[0]+1
    Fb=Psib(P[0],P[1],p,q,eta2);
    ftemp0[p]+=Fb*uh0[a];
    ftemp1[p]+=Fb*uh1[a];    
    // AB
    a=3;
    q=0;
    for(p=1;p<P[0];p++){
      Fb=Psib(P[0],P[1],p,q,eta2);
      ftemp0[p]+=Fb*uh0[a];
      ftemp1[p]+=Fb*uh1[a];
      a++;
    }
    // BC
    p=P[0];
    for(q=1;q<P[1];q++){
      Fb=Psib(P[0],P[1],p,q,eta2);
      ftemp0[p]+=Fb*uh0[a];
      ftemp1[p]+=Fb*uh1[a];
      a++;
    }
    // AC
    p=0;
    for(q=1;q<P[1];q++){
      Fb=Psib(P[0],P[1],p,q,eta2);
      ftemp0[p]+=Fb*uh0[a];
      ftemp1[p]+=Fb*uh1[a];
      a++;
    }
    // Interior modes : j runs fastest
    for(p=1;p<P[0];p++){
      for(q=1;q<P[1]-p;q++){
        Fb=Psib(P[0],P[1],p,q,eta2);
        ftemp0[p]+=Fb*uh0[a];
        ftemp1[p]+=Fb*uh1[a];
        a++;
      }
    }
    // ************************************************************************
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      aux0=0.0;
      aux1=0.0;
      for(p=0;p<Pdim;p++){
        Fa=Psia(P[0],p,eta1);
        aux0+=Fa*ftemp0[p];
        aux1+=Fa*ftemp1[p];
      }
      f0[n]=aux0;
      f1[n++]=aux1;
    }
  }
}
// ****************************************************************************
// Calcular os valores da funcao nos nos de Gauss (overloaded)
// NF = numero de campos (default = 1)
// nvar= indice da variavel (default=0)
// ****************************************************************************
void Triangle::evalGQ(double f0[],const double u0[],const int NF,const int nvar)
{
  double aux0,Fa,Fb,eta1,eta2;
  int a;
  int i,j,p,q,k;
  int n=0;
  int Pdim=P[0]+2;
  double ftemp0[Pdim];
  //printf("Triangle::evalGQ    NF = %d    nvar= %d\n",NF,nvar); 
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    // ************************************************************************
    // Construcao dos vetores temporarios
    // ************************************************************************
    for(k=0;k<Pdim;k++){
      ftemp0[k]=0.0;
    }
    // A
    a=0; p=0; q=0;
    Fb=Psib(P[0],P[1],p,q,eta2);
    ftemp0[p]+=Fb*u0[a*NF+nvar];
    // B
    a=1; p=P[0]; q=0;
    Fb=Psib(P[0],P[1],p,q,eta2);
    ftemp0[p]+=Fb*u0[a*NF+nvar];
    // C Ponto colapsado
    a=2; p=P[0]+1; q=P[1];
    Fb=Psib(P[0],P[1],p,q,eta2);
    ftemp0[p]+=Fb*u0[a*NF+nvar];   
    // AB
    a=3; q=0;
    for(p=1;p<P[0];p++){
      Fb=Psib(P[0],P[1],p,q,eta2);
      ftemp0[p]+=Fb*u0[a*NF+nvar];
      a++;
    }
    // BC
    p=P[0];
    for(q=1;q<P[1];q++){
      Fb=Psib(P[0],P[1],p,q,eta2);
      ftemp0[p]+=Fb*u0[a*NF+nvar];
      a++;
    }
    // AC
    p=0;
    for(q=1;q<P[1];q++){
      Fb=Psib(P[0],P[1],p,q,eta2);
      ftemp0[p]+=Fb*u0[a*NF+nvar];
      a++;
    }
    // Interior modes : j runs fastest
    for(p=1;p<P[0];p++){
      for(q=1;q<P[1]-p;q++){
        Fb=Psib(P[0],P[1],p,q,eta2);
        ftemp0[p]+=Fb*u0[a*NF+nvar];
        a++;
      }
    }
    // ************************************************************************
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      aux0=0.0;
      for(p=0;p<Pdim;p++){
        Fa=Psia(P[0],p,eta1);
        aux0+=Fa*ftemp0[p];
      }
      f0[n++]=aux0;
    }
  }
};
// ****************************************************************************
// Evaluates the value of the field at the vertices of the element
// ****************************************************************************
void Triangle::computeVertice(double f_vert[],const double u[], 
		    const Vertice vert[], const int map[])
{
  double aux, Fa,Fb,eta1,eta2;//x1,x2,x3;
 // double eaux,haux,aux1,faux2,faux3,gaux1,gaux2,gaux3,
 // double xa,ya,za,xb,yb,zb;//xc,yc,zc;
  int m1;
  int i,p,q;
  i=0;
  for(eta2=-1.0;eta2<=1.0;eta2+=2.0){
  
    for(eta1=eta2;eta1<=1.0;eta1+=2.0){
    
      aux=0.0;
      for(m1=0;m1<nn;m1++){
        p= mode_[m1].p_val();
        q= mode_[m1].q_val();
        Fa=Psia(P[0],p,eta1);
        Fb=Psib(P[0],P[1],p,q,eta2);
        aux+=Fa*Fb*u[m1];
      }

      f_vert[map[i++]]=aux;
    }
  }
};
// ****************************************************************************
// Evaluates the value of the field at points
// ****************************************************************************
void Triangle::computeAtPoints(const int npoints, const double LocCoord[],const double u[], 
			       const Vertice vert[], const int map[],double f[],double GloCoord[])
{
  double aux, Fa,Fb,eta1,eta2,x1,x2,x3;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  int m1;
  int i,p,q;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
  for(i=0;i<npoints;i++){
    eta1=LocCoord[2*i ];
    eta2=LocCoord[2*i+1];
 
    eaux=(1.0-eta2)/4.0;
    haux=(1.0+eta2)/2.0;
    // coordenadas (conforme Karniadakis pagina 108)
    x1=eaux*(faux1-eta1*gaux1)+haux*xc;
    x2=eaux*(faux2-eta1*gaux2)+haux*yc;
    x3=eaux*(faux3-eta1*gaux3)+haux*zc;
    GloCoord[3*i  ]=x1;
    GloCoord[3*i+1]=x2;
    GloCoord[3*i+2]=x3;
    aux=0.0;
    for(m1=0;m1<nn;m1++){
      p= mode_[m1].p_val();
      q= mode_[m1].q_val();
      Fa=Psia(P[0],p,eta1);
      Fb=Psib(P[0],P[1],p,q,eta2);
      aux+=Fa*Fb*u[m1];
    }
    //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
    f[i]=aux;
  }
};
/*
void Triangle::make_Phi(const int m,double Phi[])
{
  double eta1,eta2;
  int i,j;
  int p,q;
  int n;
  int sizeQ=Q[0]*Q[1];	
  double Fa[sizeQ],Fb[sizeQ];
  p= mode_[m].p_val();
  q= mode_[m].q_val();
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    Fb[j]=Psib(P[0],P[1],p,q,eta2);
  } 
  for(i=0;i<Q[0];i++){
    eta1=xGQ[0][i];
    Fa[i]=Psia(P[0],p,eta1);
  }
  n=0;
  for(j=0;j<Q[1];j++){
    for(i=0;i<Q[0];i++){// i runs fastest
      Phi[n++]=Fa[i]*Fb[j];
    }
  }
};
*/
void Triangle::eval_Phi(const int n,double Phi[])
{
  int q2=Q[0]*Q[1];
  for(int j=0;j<q2;j++)Phi[j]=Phi_val[n][j];
};
void Triangle::eval_GradPhi(const Vertice vert[], const int map[],const int m,double ** der)
{
  printf("Nao implementado\n");
};
// ****************************************************************************
void Triangle::Jacobian(const Vertice vert[],const int map[],double * JV)
{
  // Esta rotina eh valida so para triangulos!
  double a1,a2,a3,f1,f2,f3,g1,g2,g3;
  double xa,ya,xb,yb,xc,yc,za,zb,zc;
  double Jacobian;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;

  f1=xa-xb;
  f2=ya-yb;
  f3=za-zb;
  g1=xa-xc;
  g2=ya-yc;
  g3=za-zc;
  a1=f2*g3-f3*g2;
  a2=f3*g1-f1*g3;
  a3=f1*g2-f2*g1;
  Jacobian=sqrt(a1*a1 + a2*a2 + a3*a3)/4.0;
  if(Jacobian==0.0){
    printf("Jacobian=0\nVertices\n(%lf, %lf, %lf ), (%lf, %lf, %lf), (%lf, %lf, %lf)\n",xa,ya,za,xb,yb,zb,xc,yc,zc);
  }
//  #ifdef PRINTF_ON
//    printf("Triangle::Jacobian=%lf\n",Jacobian);
//  #endif
    for(int j=0;j<Q[1];j++)
      for(int i=0;i<Q[0];i++)
	JV[i+Q[0]*j]=Jacobian;     
};
// ****************************************************************************
/*
void Triangle::Processar_geometria(int nel,
				   const Vertice * vert,
				   const int numv, 
				   const int * VN,
				   int map[], // retorna gbnmap do elemento
				   int sgn[],
				   int sinal[],
				   int& NG,
				   int& NL,
				   std::vector<EDGE>& border,
				   int Ng[],
				   int & NF, 
				   int Face[], 
				   int Fng[], 
				   int f_mask[], 
				   const int N_add)
{
  int n0,n1,n2,sgn0,sgn1,sgn2,ng0,ng1,ng2;
//  int num,type;
  if(numv != 3){
    printf("Incompatibilidade: Dado do elemento nao eh de Triangulo\n");
    exit(0);
  }
  // passa o numero dos vertices
  n0=VN[0];
  n1=VN[1];
  n2=VN[2];
  //
  //    2
  //   / \
  //  /   \
  // 0-----1
 
  teste_aresta(P[0],n0,n1,sgn0,ng0,NG,NL,border,Ng,nel,0,vert,-1);
  teste_aresta(P[1],n1,n2,sgn1,ng1,NG,NL,border,Ng,nel,1,vert,-1);
  teste_aresta(P[1],n0,n2,sgn2,ng2,NG,NL,border,Ng,nel,2,vert, 1);
  sinal[0]=sgn0;
  sinal[1]=sgn1;
  sinal[2]=sgn2;
  //  Testar; nao criar gbnmap agora
  //n0+=N_add;
  //n1+=N_add;
  //n2+=N_add;
  //make_gbnmap(n0,n1,n2,ng0,ng1,ng2,sgn0,sgn1,sgn2,map,sgn);
 
}
*/
// ****************************************************************************

// ****************************************************************************
/*    OBSOLETA   TESTAR sem 
void Triangle::make_gbnmap(int n0,int n1,int n2,int ng0,int ng1,int ng2,
			   int sign0,int sign1,int sign2,int map[],int sgn[] )
{
  int p=P[0];
  int i;
  int imax,imin;
  int temp;
  
  //printf("SINAIS: %d %d %d\n", sign0,sign1,sign2);
  map[0]=n0; 
  sgn[0]=1;
  map[1]=n1; 
  sgn[1]=1;
  map[2]=n2; 
  sgn[2]=1;
  //printf("SINAIS de 0 1 2: %d %d %d\n", sgn[0],sgn[1],sgn[2]);
  
  imin=3;
  imax=p+1;
  // if(sign0==-1){
  temp=1;
  for(i=imin; i<=imax;i++){
    map[i]=ng0;
    sgn[i]=temp;
    temp*=sign0;// muda o sinal
    ng0++;
  }

  imin=p+2;
  imax=2*p;
 
  temp=1;
  for(i=imin;i<=imax;i++){
    map[i]=ng1;
    sgn[i]=temp;
    temp*=sign1;// muda o sinal
    ng1++;
  }
  
  imin=imax+1;
  imax=3*p-1;
  temp=1;
  for(i=imin;i<=imax;i++){
    map[i]=ng2;
    sgn[i]=temp;
    temp*=sign2;// muda o sinal
    ng2++;
  }
 
};
*/
/*
void Triangle::set_b(Mat2<double>& b, const Vertice vert[], const int map[])
{
  b.set(ndim,ndim);
  double xa,ya,za,xb,yb,zb,xc,yc,zc,eta1,eta2;
  // double x3=0.0;
  // double e3p, e3m,e2p,e2m,e1p,e1m;
  int i,j,l,m;
  double df[MAXQ*MAXQ][2];
  double a11,a12,a21,a22, J2D;
  double aux,aux0,aux1;
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  // calculo dos coeficientes aij
  a11=0.5*(xb-xa);
  a12=0.5*(xc-xa);
  a21=0.5*(yb-ya);
  a22=0.5*(yc-ya);
  // calculo de J2D
  J2D=a11*a22 - a12*a21;
  // calculo da matriz b[i][j]
  b[0][0]=a22/J2D;
  b[0][1]=-a12/J2D;
  b[1][0]=-a21/J2D;
  b[1][1]=a11/J2D;
};
*/

// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// ****************************************************************************
// Necessita ser reescrito para considerar o caso do triangulo no espaco 3d

void Triangle::Gradiente(double * grad[],
												 const  double fvec[],
												 const Vertice vert[], const int map[])
{
  double xa,ya,xb,yb,xc,yc,eta1,eta2;
 // double x3=0.0;
 // double e3p, e3m,e2p,e2m,e1p,e1m;
  int i,j,l,m;
  double df[MAXQ*MAXQ][2];
  double a11,a12,a21,a22, J2D;
  double b[2][2];
  double aux,aux0,aux1;
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
 // za=vert[map[0]].z;
 // zb=vert[map[1]].z;
 // zc=vert[map[2]].z;
  // calculo dos coeficientes aij
  a11=0.5*(xb-xa);
  a12=0.5*(xc-xa);
  a21=0.5*(yb-ya);
  a22=0.5*(yc-ya);
  // calculo de J2D
  J2D=a11*a22 - a12*a21;
  // calculo da matriz b[i][j]
  b[0][0]=a22/J2D;
  b[0][1]=-a12/J2D;
  b[1][0]=-a21/J2D;
  b[1][1]=a11/J2D;
  
  // calculo das derivadas com relacao a eta1 e eta2
  for(int q=0;q<Q[1];q++){
    for(int p=0;p<Q[0];p++){
      m=p+q*Q[0];
      aux0=0.0;
      aux1=0.0;
      for(l=0;l<Q[0];l++){
        aux0+=D[p][l][0]*fvec[(l+q*Q[0])];
        aux1+=D[q][l][1]*fvec[(p+l*Q[0])];
      }
      df[m][0]=aux0;
      df[m][1]=aux1;
    }
  } // Q[1]*Q[0]*Q[0]  operacoes
  
  // calculo das derivadas com relacao a xi1, xi2
  for(int q=0;q<Q[1];q++){
    eta2=xGQ[1][q];
    for(int p=0;p<Q[0];p++){
      eta1=xGQ[0][p];
      m=p+q*Q[0];
      aux=(1-eta2);
      aux0=(2/aux)*df[m][0]; // d/dxi_1
      aux1=(1+eta1)/aux*df[m][0] + df[m][1];// d/dxi_2
      df[m][0]=aux0;// d/dxi_1
      df[m][1]=aux1;// d/dxi_2
    }
  } // Q[1]*Q[0]  operacoes
  
  // calculo das derivadas com relacao a x1 e x2
  for(int q=0;q<Q[1];q++){
    for(int p=0;p<Q[0];p++){
      m=p+q*Q[0];
      for(i=0;i<2;i++){
        aux=0;
        for(j=0;j<2;j++){
          aux+=b[j][i]*df[m][j];
        }
        grad[i][m]=aux;
      }
    }
  } // 4*Q[1]*Q[0]  operacoes
};
// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// Checa o gradiente e imprime o teste no arquivo fout
// ****************************************************************************
// Necessita ser reescrito para considerar o caso do triangulo no espaco 3d

void Triangle::Gradiente(FILE * fout, double * grad[],
												 const  double fvec[],
												 const Vertice vert[], const int map[])
{
  double xa,ya,xb,yb,xc,yc,eta1,eta2,x1,x2;//za,zb,zc;
  double x3=0.0;
  double e2p,e2m,e1p,e1m;
  int i,j,l,m;
  double df[MAXQ*MAXQ][2]; // 2 ou 3 ?
  double a11,a12,a21,a22, J2D;
  double b[2][2];
  double aux,aux0,aux1;
 
	// coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
 // za=vert[map[0]].z;
 // zb=vert[map[1]].z;
 // zc=vert[map[2]].z;
  // calculo dos coeficientes aij
  a11=0.5*(xb-xa);
  a12=0.5*(xc-xa);
  a21=0.5*(yb-ya);
  a22=0.5*(yc-ya);
  // calculo de J2D
  J2D=a11*a22 - a12*a21;
  // calculo da matriz b[i][j]
  b[0][0]=a22/J2D;
  b[0][1]=-a12/J2D;
  b[1][0]=-a21/J2D;
  b[1][1]=a11/J2D;
  
  // calculo das derivadas com relacao a eta1 e eta2
  for(int q=0;q<Q[1];++q){
    for(int p=0;p<Q[0];++p){
      m=p+q*Q[0];
      aux0=0.0;
      aux1=0.0;
      for(l=0;l<Q[0];++l){
        aux0+=D[p][l][0]*fvec[(l+q*Q[0])];
        aux1+=D[q][l][1]*fvec[(p+l*Q[0])];
      }
      df[m][0]=aux0;
      df[m][1]=aux1;
    }
  }
  
  // calculo das derivadas com relacao a xi1, xi2
  for(int q=0;q<Q[1];++q){
    eta2=xGQ[1][q];
    for(int p=0;p<Q[0];++p){
      eta1=xGQ[0][p];
      m=p+q*Q[0];
      aux=(1-eta2);
      aux0=(2/aux)*df[m][0]; // d/dxi_1
      aux1=(1+eta1)/aux*df[m][0] + df[m][1];// d/dxi_2
      df[m][0]=aux0;// d/dxi_1
      df[m][1]=aux1;// d/dxi_2
    }
  }
  // calculo das derivadas com relacao a x1 e x2
  for(int q=0;q<Q[1];q++){
    for(int p=0;p<Q[0];p++){
      m=p+q*Q[0];
      for(i=0;i<2;i++){
        aux=0;
        for(j=0;j<2;j++){
          aux+=b[j][i]*df[m][j];
        }
        grad[i][m]=aux;
      }
    }
  }
  
  // Cheque do gradiente
  for(j=0;j<Q[1];++j){
    eta2=xGQ[1][j];
    e2p=(1.0+eta2)/2.0;
    e2m=(1.0-eta2)/2.0;
    for(i=0;i<Q[0];++i){
      eta1=xGQ[0][i];
      e1p=(1.0+eta1)/2.0;
      e1m=(1.0-eta1)/2.0;
      // coordenadas x1, x2, x3
      x1=(xa*e1m+xb*e1p)*e2m+e2p*xc;
      x2=(ya*e1m+yb*e1p)*e2m+e2p*yc;
      
      m=i+j*Q[0];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,grad[0][m],grad[1][m],g(x1,x2,x3));
    }
  }
  
};

// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// ****************************************************************************

void Triangle::Gradiente(FILE * fout, double * grad[],
												 double (*func)(double, double, double),
												 const Vertice vert[], const int map[])
{
  double xa,ya,xb,yb,xc,yc,eta1,eta2,x1,x2;//za,zb,zc;
  double x3=0.0;// triangulo no plano x,y
  double e2p,e2m,e1p,e1m;
  int i,j,m;
  double fvec[MAXQ*MAXQ];
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
 // za=vert[map[0]].z;
 // zb=vert[map[1]].z;
 // zc=vert[map[2]].z;
  
  // calculo do vetor contendo a funcao nos pontos de integracao de Gauss
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    e2p=(1.0+eta2)/2.0;
    e2m=(1.0-eta2)/2.0;
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      e1p=(1.0+eta1)/2.0;
      e1m=(1.0-eta1)/2.0;
      // coordenadas x1, x2, x3
      x1=(xa*e1m+xb*e1p)*e2m+e2p*xc;
      x2=(ya*e1m+yb*e1p)*e2m+e2p*yc;
      
      m=i+j*Q[0];
   
      fvec[m]=func(x1,x2,x3);
    }
  }
  
  Gradiente(fout,grad,fvec,vert,map);
  
};
// ****************************************************************************
void Triangle::print_nome(FILE * fout)
{
  fprintf(fout,"ELEMENTO TRIANGULAR\n");
};
// ****************************************************************************
// Imposes Dirichlet type boundary conditions
// ****************************************************************************
void Triangle::Dirichlet(const int aresta,
                         const Vertice vert[],
                         const int vert_map[],
                         const int nmap[],
                         const int sgn[],
                         int bflag[],
                         double X[],
                         double (*f)(double,double,double))
{
  //cout << "Entrando em Triangle::Dirichlet" << endl;
  // **************************************************************************
  // flag = 0 : Dirichlet, valor conhecido, bflag=0
  //      = 1 : valor desconhecido, bflag=1
  // **************************************************************************
  int flag=0;
  int i,j,k,p,q,nd,ii;
  double xa,ya,za,xb,yb,zb;
  double eta1,x1,x2,x3;
  double gaux1,gaux2,gaux3;
  
  double aux,a0,ap;
  double J=0.0;  
  
  if(aresta==0)//aresta 0; nos 0 e 1
    {	
      xa=vert[vert_map[0]].x;
      xb=vert[vert_map[1]].x;
      ya=vert[vert_map[0]].y;
      yb=vert[vert_map[1]].y;
      za=vert[vert_map[0]].z;
      zb=vert[vert_map[1]].z;
      nd=0;//direcao 0
    }
  else if(aresta==1)// aresta 1; nos 1 e 2
    {
      xa=vert[vert_map[1]].x;
      xb=vert[vert_map[2]].x;
      ya=vert[vert_map[1]].y;
      yb=vert[vert_map[2]].y;
      za=vert[vert_map[1]].z;
      zb=vert[vert_map[2]].z;
      nd=1;//direcao 1
    }
  else// aresta 2; nos 0 e 2
    {
      xa=vert[vert_map[0]].x;
      xb=vert[vert_map[2]].x;
      ya=vert[vert_map[0]].y;
      yb=vert[vert_map[2]].y;
      za=vert[vert_map[0]].z;
      zb=vert[vert_map[2]].z;
      nd=1;//direcao 1
    }
  p=P[nd];
  q=Q[0];
  double x[q],wtemp[q],Dtemp[MAXQ][MAXQ];
  Gauss_Lobatto_Jacobi_parameters(q,0.0,0.0,x,wtemp,Dtemp);
  double func[q];
  double psi[p+1][q];
  int p1=p-1;
 
//  int Ti[p1*p1],Tj[p1*p1];
//  double Tx[p1*p1];
//  int count=0;

#ifdef _NEWMAT
  NEWMAT::Matrix A(p1,p1);
  NEWMAT::ColumnVector B(p1), Y(p1);
// #else 
//   double A[p1][p1];
//   double B[p1],Y[p1];
#endif

  a0=f(xa,ya,za);
  ap=f(xb,yb,zb);

  if(p>1){// p>1
    gaux1=xa-xb;
    gaux2=ya-yb;
    gaux3=za-zb;
    J=sqrt(gaux1*gaux1+gaux2*gaux2+gaux3*gaux3)/2.0;
    
    for(k=0;k<q;k++)
      {
        eta1= x[k];
        aux=(eta1+1.0)*0.5;
        x1=xa-aux*gaux1;
        x2=ya-aux*gaux2;
        x3=za-aux*gaux3;
        for(i=0;i<=p;i++)
        {
          psi[i][k]=Psia(p,i,eta1);
        }
          func[k]=f(x1,x2,x3);
      }
      for(i=1;i<p;i++){
        for(j=i;j<p;j++){
          aux=0.0;
          for(k=0;k<q;k++) aux+=psi[i][k]*psi[j][k]*wtemp[k];
          aux*=J;
// 	Ti[count]=i-1;
// 	Tj[count]=j-1;
// 	Tx[count++]=aux;
          A.element(i-1,j-1)=aux;
          if(j!=i){
// 	  Ti[count]=j-1;
// 	  Tj[count]=i-1;
// 	  Tx[count++]=aux;
            A.element(j-1,i-1)=aux;
          }
        }
        aux=0.0;
        for(k=0;k<q;k++)
          aux+=psi[i][k]*(func[k]-a0*psi[0][k]-ap*psi[p][k])*wtemp[k];
          B.element(i-1)=aux*J;
      }
    
#ifdef _NEWMAT
    Y = A.i() * B; B=Y;
// #else 
//     ResolveSistema(p1,count,Ti,Tj,Tx,B,Y);
#endif
  }
  // mapeamento do resultado no vector global
  int temp;
  if(aresta==0){
    //temp=gbnmap[0]*NFields+varn;
    temp=nmap[0];
    X[temp]=sgn[0]*a0;
    bflag[temp]=flag;
    //temp=gbnmap[1]*NFields+varn;
    temp=nmap[1];
    X[temp]=sgn[1]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++){ 
      ii=i+2;
      //temp=gbnmap[ii]*NFields+varn;
      temp=nmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  else if(aresta==1){
    //temp=gbnmap[1]*NFields+varn;
    temp=nmap[1];
    X[temp]=sgn[1]*a0;
    bflag[temp]=flag;
    //temp=gbnmap[2]*NFields+varn;
    temp=nmap[2];
    X[temp]=sgn[2]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++) {
      ii=i+P[0]+1;
      //temp=gbnmap[ii]*NFields+varn;
      temp=nmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  else{
    //temp=gbnmap[0]*NFields+varn;
    temp=nmap[0];
    X[temp]=sgn[0]*a0;
    bflag[temp]=flag;
    //temp=gbnmap[2]*NFields+varn;
    temp=nmap[2];
    X[temp]=sgn[2]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++) {
      ii=i+P[0]+P[1];
      //temp=gbnmap[ii]*NFields+varn;
      temp=nmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  //cout << "Saindo de Triangle::Dirichlet" << endl;
};

void Triangle::teste(int & v)
{
   v=100000;
};


// ****************************************************************************
// Evaluates the value of the field at the Gauss Quadrature points
// ****************************************************************************
void Triangle::computeFuncGQ(double f_[], 
			     const Vertice vert[], const int map[],
			     double (*func)(double,double,double))
{
  double eta1,eta2,x1,x2,x3;
  double eaux,faux1,faux2,faux3,gaux1,gaux2,gaux3,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc;
  int i,j,count;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  faux1=xa+xb;
  faux2=ya+yb;
  faux3=za+zb;
  gaux1=xa-xb;
  gaux2=ya-yb;
  gaux3=za-zb;
  count=0;
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/4.0;
    haux=(1.0+eta2)/2.0;
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      // coordenadas (conforme Karniadakis pagina 108)
      x1=eaux*(faux1-eta1*gaux1)+haux*xc;
      x2=eaux*(faux2-eta1*gaux2)+haux*yc;
      x3=eaux*(faux3-eta1*gaux3)+haux*zc;
      f_[count++]=func(x1,x2,x3);
    }
  }
};
// ****************************************************************************
//   void Triangle::localFaceModeMap(const int fnum,
//   			  const int gbnmap[], int bflag[],
//   			  int lmap[],
//   			  const int varn, const int NFields)
//   {
//     int temp;
//     int a,i,j,nmax;
//     double b[MAXNB];
//     int flag=0;
//     if(fnum==0){ //AB
//       bflag[gbnmap[0] * NFields+varn ]=flag;// A
//       bflag[gbnmap[1] * NFields+varn ]=flag;// B
//       a=3;
//       j=0;
//       nmax=P[0];
//       lmap[0]=0;
//       lmap[nmax]=1;
//       for(i=1;i<nmax;i++){
//         bflag[gbnmap[a] * NFields+varn ]=flag;
//         lmap[i]=a;
//         a++;
//       }
//     }
//     else if(fnum==1){ //BC 
//       bflag[gbnmap[1] * NFields+varn ]=flag;// B
//       bflag[gbnmap[2] * NFields+varn ]=flag;// C
//       a=P[0]+2;
//       i=P[0];
//       nmax=P[1];
//       lmap[0]=1;
//       lmap[nmax]=2;
//       for(int j=1;j<nmax;j++){
//         bflag[gbnmap[a] * NFields+varn ]=flag;
//         lmap[j]=a;
//         a++;
//       }
//     }
//     else if(fnum==2){ //AC 
//       bflag[gbnmap[0] * NFields+varn ]=flag;// A
//       bflag[gbnmap[2] * NFields+varn ]=flag;// C
//       a=2*P[0]+1;
//       i=0;
//       nmax=P[1];
//       lmap[0]=0;
//       lmap[nmax]=2;
//       for(int j=1;j<nmax;j++){
//         bflag[gbnmap[a] * NFields+varn ]=flag;
//         lmap[j]=a;
//         a++;
//       }
//     }
//   };
/*
void Triangle::printStdel() const
{
  printf("ndim= %d\nP   Q  gqt\n",ndim);
  for (int i=0; i<ndim; i++)
    printf("%d   %d   %d\n",P[i], Q[i], gqt[i]);
};

*/
// ****************************************************************************
/*
void Triangle::make_mass(double ** MM, const double JV[])
{
  int Nb,ni; 
  int i,ii,j,jj;
  ni=nn-nb;
#ifdef PRINTF_ON
  printf("Making local matrices 1: nb = %d  Ni= %d\n", nb,Ni);
#endif
  int p,q,r,s;
  int flag1, flag2;
  double aux, aaux, Fa,Fb,Ga,Gb;
  int m1,m2;
  for(m1=0;m1<nn;m1++){
    // indexes of point m1
    p= mode_[m1].p_val();
    q= mode_[m1].q_val();
    for(m2=m1;m2<nn;m2++){
      // indexes of point m2
      r= mode_[m2].p_val();
      s= mode_[m2].q_val();
      aux=0.0;
      for(j=0;j<Q[1];j++){
	aaux=0.0;
	for(i=0;i<Q[0];i++){
	  Fa=Psia(P[0], p,xGQ[0][i]);
	  Ga=Psia(P[0], r,xGQ[0][i]);
	  aaux+=wGQ[0][i]*Fa*Ga*JV[i+j*Q[0]];//<===Jacobiano======
	}
	Fb=Psib(P[0], P[1], p, q, xGQ[1][j]);
	Gb=Psib(P[0], P[1], r, s, xGQ[1][j]);
	aux+=wGQ[1][j]*Fb*Gb*aaux;
      }
      MM[m1][m2]=aux;
      MM[m2][m1]=aux;
    }
  }
};
*/
//   void Triangle::Jacobian_Vector(const Vertice vert[],const int map[],double JV[])
//   {
//     int q2=Q[0]*Q[1];
//     double J=Jacobian(vert,map);
//     for(int i=0;i<q2;i++)JV[i]=J;
//   };


// *****************************************************************************
// Calcula os tracos de Phi, GradPhi e Jb nos pontos de Gauss sobre as arestas
// ordenando-os de forma que haja coincidencia de pontos dos elementos vizinhos
// *****************************************************************************
void Triangle::elem_traces(const Vertice vert[],const int map[],const int sinal[],
			   double *** TP,double **** TGP,double * Jb)
{
  double xa,ya,xb,yb,xc,yc,eta1,eta2;
 // double x3=0.0;
  int h,l,m;
  double a11,a12,a21,a22, J2D;
  double b[2][2];
  double d1,d2,aux0,aux1,aux2,der1,der2;  
// coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  // calculo dos coeficientes aij
  a11=0.5*(xb-xa);
  a12=0.5*(xc-xa);
  a21=0.5*(yb-ya);
  a22=0.5*(yc-ya);
  // calculo de J2D
  J2D=a11*a22 - a12*a21;
  // calculo da matriz b[i][j]
  b[0][0]=a22/J2D;
  b[0][1]=-a12/J2D;
  b[1][0]=-a21/J2D;
  b[1][1]=a11/J2D;
  
  //pontos de Gauss em uma dimensao
  double x[qborder], wtemp[qborder], Dtemp[MAXQ][MAXQ];
  Gauss_Jacobi_parameters(qborder, 0.0, 0.0, x, wtemp, Dtemp);
  // Gauss_Jacobi integra precisamente P_(2*q-1)
  // Obs.: Gauss_Lobatto_Jacobi da problemas pois calcula o valor
  // de Phi no ponto colapsado. Alem disso,
  // Gauss_Lobatto_Jacobi integra precisamente P_(2*q-3): 
  // grau menor que Gauss_Jacobi
  h=0;
  aux0=sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya)) / 2.0;
  for(int l=0;l<qborder;l++){
    Jb[h*qborder+l]=aux0*wtemp[l];
  }
  h=1;
  aux0=sqrt( (xc-xb)*(xc-xb) + (yc-yb)*(yc-yb)) / 2.0;
  for(int l=0;l<qborder;l++){
    Jb[h*qborder+l]=aux0*wtemp[l];
  }
  h=2;
  aux0=sqrt( (xc-xa)*(xc-xa) + (yc-ya)*(yc-ya)) / 2.0;
  for(int l=0;l<qborder;l++){
    Jb[h*qborder+l]=aux0*wtemp[l];
  }
 
 // int s0=nn*ndim*qborder;
 // int s1=   ndim*qborder;
  
  for(m=0;m<nn;m++) { // loop sobre os modos
    int p= mode_[m].p_val();
    int q= mode_[m].q_val();
    
    // aresta  h = 0
    h=0;
   //int i_m = m*qborder;
    eta2=-1.0;
    aux2=Psib(P[0],P[1],p,q,eta2);
    d2=DPsib(P[0],P[1],p,q,eta2);
    for(int l=0;l<qborder;l++){ // loop sobre os pontos de Gauss
      eta1=sinal[h]*x[l];
      aux1=Psia(P[0],p,eta1);
      d1=DPsia(P[0],p,eta1);
      /*TP[i_m + l]*/ TP[h][m][l]=aux1*aux2;
      der1 = d1*aux2;
      der2 = 0.5*(1.0+eta1)*der1 + aux1*d2;
      /*TGP[h*s0+m*s1     +l]*/    TGP[h][m][0][l] = der1 * b[0][0] + der2 * b[1][0];
      /*TGP[h*s0+m*s1+qborder+l]*/ TGP[h][m][1][l] = der1 * b[0][1] + der2 * b[1][1];
    } // loop sobre os pontos de Gauss
    
    // aresta  h = 1
    h=1;
  //  i_m = nn*qborder + m*qborder;
    eta1=1.0;
    aux1=Psia(P[0],p,eta1);
    d1=DPsia(P[0],p,eta1);
    for(int l=0;l<qborder;l++){ // loop sobre os pontos de Gauss
      eta2=sinal[h]*x[l];
      aux2=Psib(P[0],P[1],p,q,eta2);
      d2=DPsib(P[0],P[1],p,q,eta2);
      /*TP[i_m + l]*/ TP[h][m][l]=aux1*aux2;
      //TP[i_m + l]=aux1*aux2;
      der1 = 2.0/(1-eta2)*d1*aux2;
      der2 = der1 + aux1*d2;
      /*TGP[h*s0+m*s1     +l]*/    TGP[h][m][0][l] = der1 * b[0][0] + der2 * b[1][0];
      /*TGP[h*s0+m*s1+qborder+l]*/ TGP[h][m][1][l] = der1 * b[0][1] + der2 * b[1][1];
    } // loop sobre os pontos de Gauss
    
    // aresta  h = 2
    h=2;
   // i_m = h*nn*qborder + m*qborder;
    eta1=-1.0;
    aux1=Psia(P[0],p,eta1);
    d1=DPsia(P[0],p,eta1);
    for(l=0;l<qborder;l++){ // loop sobre os pontos de Gauss
      eta2=sinal[h]*x[l];
      aux2=Psib(P[0],P[1],p,q,eta2);
      d2=DPsib(P[0],P[1],p,q,eta2);
      /*TP[i_m + l]*/ TP[h][m][l]=aux1*aux2;
      //TP[i_m + l]=aux1*aux2;
      der1 = 2.0/(1-eta2)*d1*aux2;
      der2 = aux1*d2;
      /*TGP[h*s0+m*s1     +l]*/    TGP[h][m][0][l] = der1 * b[0][0] + der2 * b[1][0];
      /*TGP[h*s0+m*s1+qborder+l]*/ TGP[h][m][1][l] = der1 * b[0][1] + der2 * b[1][1];
    } // loop sobre os pontos de Gauss
    
  } // loop sobre os modos
  //for(int i=0;i<nborder*nn*ndim*qborder;i++)printf("Triangle TGP[%d]= %g\n",i,TGP[i]);
}

//#undef PRINTF_ON
// ****************************************************************************
// Calcula o traco e o coloca na ordem correta de acordo com o sinal da borda *
// ****************************************************************************
void Triangle::trace(const int lado, const int qmax, const int sinal, 
		     const double * valores, // valores nos pontos de Gauss
		     double * saida)
{
  //int qmax = qborder;
  int nd,ind,inc;
  
  if(lado == 0){
		nd=0;
    ind=0;
    inc=1;
  }
  else {
		nd=1;
    if (lado == 1) {
      ind=Q[0]-1;
      inc=Q[0];
    }
    else {
      ind=0;
      inc=Q[0];
    }
  }
	
	int q=Q[nd];
  double temp[q];
	
  for(int i = 0; i < q; ++i) {
    saida[i] = valores[ind];
    ind += inc;
  }
	// Versao anterior era assim
  // ******************************************************************
  // Fazer interpolacao de Lagrange sobre novos pontos de integracao  *
  // Quando for a direcao 1 (nd=1) onde ha o vertice colapsado ou     *
  // quando o numero de pontos de quadratura q for menor que o numero *
  // de elementos do vetor de traco. Exclusivo para Triangle !!!!     *
  // ******************************************************************

    double Old[q];

    for(int i = 0; i < q; ++i) {
      temp[i]=saida[i];
      Old[i]=xGQ[nd][i];
    }

    for(int k=0;k<q;++k){
      double prod=1.0;
      double y=Old[k];
      for(int j=0;j<q;++j){
				if(j!=k)prod*=(y-Old[j]);
      }
      temp[k]/=prod;
    }
		// ****************************************************************
		// Os pontos de Gauss onde estao os tracos sao do tipo
		// Gauss-Jacobi (nao incluimos os pontos extremos),
    // portanto diferentes dos pontos de Gauss
		// do elemento. Logo toda aresta necessita calcular o traco
		// ****************************************************************
    double Jac[qmax],wtemp[qmax],Dtemp[MAXQ][MAXQ];
    Gauss_Jacobi_parameters(qmax,0.0,0.0,Jac,wtemp,Dtemp);

    for(int i=0; i < qmax; ++i) {
      double sum=0.0;
      double y=Jac[i] * sinal; // percorre os pontos de Gauss na ordem decrescente
      for(int k=0;k<q;++k) {   // se o sinal == -1
				double prod=1.0;
				for(int j=0;j<q;++j) {
					if(j!=k)prod*=(y-Old[j]);
				}
				sum+=(temp[k]*prod);
      }
      saida[i]=sum;
    }
  
};
// revisado em 25/10/2011

const int Triangle::aresta_lvert(const int & i, const int & j) const {return aresta[i][j];};
const int Triangle::face_lvert(const int & i, const int & j) const {return face[i][j];};
const int Triangle::show_nvf(const int &i) const {return 0;};
const int Triangle::show_face_tipo(const int &i) const {return 0;};
const int Triangle::show_fd0(const int &i) const {return 0;};
const int Triangle::show_fd1(const int &i) const {return 0;};
const int Triangle::show_fv2(const int &i) const {return 0;};
const int Triangle::show_ind_mode(const int & i, const int & j, const int & k) const {return ind_mode_[i][j];};// a ser implementado
void Triangle::superficie_externa(const int Vert_map[], const Vertice vert[],
                                  const int & num_local,
                                  double & area,double normal[3])
{
  int v0,v1;
  v0=Vert_map[aresta[num_local][0]];
  v1=Vert_map[aresta[num_local][1]];
  double lx,ly,lz;
  lx=vert[v1].x - vert[v0].x;
  ly=vert[v1].y - vert[v0].y;
  lz=vert[v1].z - vert[v0].z;
  
  normal[0] = -sinal_normal[num_local] * ly;
  normal[1] =  sinal_normal[num_local] * lx;
  normal[2] = 0.0;
  
  area = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
  
  normal[0]/=area;
  normal[1]/=area;
};
