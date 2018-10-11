#include "spectral.h"
//#define PRINTF_ON
// ****************************************************************************
// Class Linear
// ****************************************************************************
// ****************************************************************************
Linear::Linear(int p0,int q0) // 10/02/2013
{
  //  printf("Iniciando Linear\n");
  tipo=1;
  ndim=1;
  nv=2;
  ne=0;
  nf=0;
  nb=2;
  nborder=2;
  qborder = 1;
  vtk_type=3;
  emapi = new int [nborder+1];
  emapv = new int [nborder];
  set(p0, q0);
};
// ****************************************************************************
Linear::~Linear()  // 10/02/2013
{
  //  printf("Destruindo Linear\n");
  delete [] emapi; emapi=nullptr;
  delete [] emapv; emapv=nullptr;
  delete [] D_Phi_val; D_Phi_val=nullptr;

  //libera memoria de ind_mode_
  delete [] ind_mode_;  ind_mode_ = nullptr;
};
// ****************************************************************************

// Inicializa os dados do elemento
void Linear::set(int p0, int q0) // 10/02/2013
{
  int i;
  // Mapeamento das arestas
  // Indices que marcam o inicio das arestas
  emapi[0]=0;
  emapi[1]=1;
  emapi[2]=2;
  emapv[0]=0;// Vertice A
  emapv[1]=1;// Vertice B

  //printf("Em Linear::set(%d,%d)\n",p0,q0);
  P[0]=p0; Q[0]=q0; gqt[0]=3;//Gauss-Lobatto-Jacobi (inclui os dois vertices)
  P[1]=1;  Q[1]=1;  gqt[1]=0;//Segunda dimensao
  P[2]=1;  Q[2]=1;  gqt[2]=0;// Terceira dimensao
  NGQP=Q[0];
  //aloca memoria para ind_mode_
  ind_mode_ = new int [P[0]+1];
#ifdef PRINTF_ON
  printf("Em Linear::set apos Sair de Linear::setStdel\n");
  printf("P: %d %d\nQ: %d %d\ngqt : %d %d\n",P[0],P[1],Q[0],Q[1],gqt[0],gqt[1]);
#endif

  // Mapeamento de nos, modos e arestas
  mode_[0].set_mode(0,0,0);ind_mode_[0]=0;// A
  mode_[1].set_mode(P[0],0,0);ind_mode_[1]=1;// B
  // AB
  int a=2;
  for(i=1;i<P[0];i++){
    mode_[a].set_mode(i,0,0);ind_mode_[i]=a;
    a++;
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
  printf("Saindo de Linear::set: nv = %d  nn = %d (< MAXMODES = %d) nb = %d (< MAXNB = %d)\n",nv,nn,MAXMODES,nb,MAXNB);
 printf("Linear::set nn = %d nb= %d\n",nn,nb);
#endif

  // Construcao da matriz Phi_val[nn][Q[0]]
  double eta1;
  for(int m=0;m<nn;m++){
    int p= mode_[m].p_val();

    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      Phi_val[m][i]=Psia(P[0],p,eta1);
    }
  }
};
// 10/02/2013
// ****************************************************************************
void Linear::print(FILE * fout)  // 10/02/2013
{
  Stdel::printStdel();
  fprintf(fout,"Linear::print\n");
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

double Linear::mass(int m1, int m2, const double JV[]) // 10/02/2013
{
  int p1,p2;
  int i;
  double aux,Fa,Ga;
  //printf("\nEntrou Linear::mass(%d, %d)\n",m1,m2);
  // indexes of point m1
  p1= mode_[m1].p_val();
  // indexes of point m2
  p2= mode_[m2].p_val();
  aux=0.0;
  for(i=0;i<Q[0];i++){
    Fa=Psia(P[0],p1,xGQ[0][i]);//primeiro
    Ga=Psia(P[0],p2,xGQ[0][i]);
    aux+=wGQ[0][i]*Fa*Ga;//*JV[i+Q[0]*j];// <==multiplica pelo Jacobiano-<
  }
  //printf("Saiu Linear::mass(%d, %d) = %g\n",m1,m2,aux*JV[0]);
  return (aux*JV[0]);// <----------multiplica pelo Jacobiano----<
}
// ****************************************************************************
// Creates the local matrices using static condensation
// It uses NEWMAT matrices and operators
// ****************************************************************************
void Linear::make_local_matrices()  // 10/02/2013
{
  int ni;
  int i,ii,j,jj;
  ni=nn-nb;
#ifdef PRINTF_ON
  printf("Making local matrices 1: nb = %d  ni= %d\n", nb,ni);
#endif

  double M[nn][nn];
  int p1,p2;
  double aux, Fa,Ga;
  int m1,m2;
  for(m1=0;m1<nn;m1++){
    // indexes of point m1
    p1= mode_[m1].p_val();
    for(m2=m1;m2<nn;m2++){
      // indexes of point m2
      p2= mode_[m2].p_val();
      aux=0.0;
      for(i=0;i<Q[0];i++){
	Fa=Psia(P[0], p1,xGQ[0][i]);
	Ga=Psia(P[0], p2,xGQ[0][i]);
	aux+=wGQ[0][i]*Fa*Ga;
      }
      M[m1][m2]=aux;
      M[m2][m1]=aux;
    }
  }
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
#ifdef PRINTF_ON
  printf("Saindo Linear::make_local_matrices: nb = %d  ni= %d\n", nb,ni);
#endif
};

// ****************************************************************************

void Linear::make_mass_matrices(int NFields)
{/*
 //printf("Linear make_mass_matrices\n");
 // ***************************
 make_local_matrices();
 // ***************************
 if(NFields > 1)duplicar_mass_matrix(NFields);
 #ifdef PRINTF_ON
 printf("Saindo Linear::make_mass_matrices");
 printf("(%2d): nb = %d  ni= %d\n",NFields,nb,nn-nb);
 #endif
 */
};

// ****************************************************************************
// Sets Gauss parameters for the triangle
// Ja em coordenadas eta1 e eta2, incluindo o termo (1-eta2)/2
// (Gauss-Jacobi quadratura) Ver Apendice B de Karniadakis & Sherwin
// ****************************************************************************
void Linear::gauss_parameters_default()  // 10/02/2013
{
  double x[MAXQ];
  double wtemp[MAXQ];
  double Dtemp[MAXQ][MAXQ];

  // Valores para a segunda e terceira dimensoes: 2 pontos (eta=-1.0) com w=1.0
  xGQ[1][0]= -1.0;
  wGQ[1][0]= 1.0;
  xGQ[2][0]= -1.0;
  wGQ[2][0]= 1.0;

  switch (gqt[0]) {

  case 1:
    Gauss_Jacobi_parameters(Q[0], 0.0, 0.0, x, wtemp, Dtemp);
    break;

  case 2:
    Gauss_Radau_Jacobi_parameters(Q[0], 0.0, 0.0, x, wtemp, Dtemp);
    break;

  case 3:
    Gauss_Lobatto_Jacobi_parameters(Q[0], 0.0, 0.0, x, wtemp, Dtemp);
    break;
  }

  // printf("Linear: Gauss parameters defaults: tipo %d\n",gqt[0]);
  for(int j=0; j<Q[0]; j++){
    xGQ[0][j] = x[j];
    wGQ[0][j] = wtemp[j];
    // printf(" Linear wGQ[0][%d] = %g\n", j,wGQ[0][j]);
  }

  for(int k=0;k<Q[0];k++) {
    for(int l=0;l<Q[0];l++){
      D[k][l][0]=Dtemp[k][l];
      //printf(" Linear D[%d][%d][0] = %g\n", k,l,D[k][l][0]);
    }
  }
};
// 10/02/2013
// ****************************************************************************
// OBSERVACAO: MULTIPLICAR POR JV
// ****************************************************************************
void Linear::vector_of_integral_of_f_Phi_dv(double vec[],
				    double (*func)(double, double, double),
				    const Vertice vert[], const int map[],
				    const double JV[])  // 10/02/2013
{
  int n,p;
  int i;
  double aux, Fa,eta1,x1;
  double eaux,aux1;
  double xa,xb;
  double x2=0.0;
  double x3=0.0;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  aux1=xb-xa;

  for(n=0;n<nn;n++){
    p= mode_[n].p_val();
    aux=0.0;
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      // coordenadas
      eaux=(eta1+1.0)/2.0;
      x1=eaux*aux1+xa;
      Fa=Psia(P[0],p,eta1);
      aux+=(wGQ[0][i]*Fa*func(x1,x2,x3));//*JV[i+Q[0]*j];// <======Jacobiano==<
    }
    // ------Multiplicar pelo Jacobiano----
    vec[n]=(aux*JV[0]);
  }
};
// ****************************************************************************
void Linear::vector_of_integral_of_f_Phi_dv(double vec[],
				    const double func[],
				  //  const Vertice vert[],const int map[],
				    const double JV[])  // 10/02/2013
{
  int n,p;
  int i;
  double aux, Fa,eta1;//x1;
//  double eaux,aux1;
//  double xa,xb;
//  xa=vert[map[0]].x;
//  xb=vert[map[1]].x;

//  aux1=xb-xa;

  for(n=0;n<nn;n++){
    p= mode_[n].p_val();
    aux=0.0;
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      // coordenadas
//      eaux=(eta1+1.0)/2.0;
    //  x1=eaux*aux1+xa;
      Fa=Psia(P[0],p,eta1);
      aux+=(wGQ[0][i]*Fa*func[i]);//*JV[i+Q[0]*j];// <===Jacobiano==<
    }
    // ------Multiplicar pelo Jacobiano----
    vec[n] = aux*JV[0];// <------Multiplicar pelo Jacobiano----<
  }
};
// ****************************************************************************

// ****************************************************************************
void Linear::printtofile(FILE * fout,const double u[],
			 double (*func)(double,double,double),
			 const Vertice vert[], const int map[])  // 10/02/2013
{
  double aux,x1,x2,x3;
  double eaux,aux1;
  double xa,xb;
  double ftemp[Q[0]];
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  x2=0.0;
  x3=0.0;

  aux1=xb-xa;

  int m1=0;
  evalGQ(ftemp,u);
  for(int i=0;i<Q[0];i++){
    // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
    eaux=(xGQ[0][i]+1.0)/2.0;
    x1=eaux*aux1+xa;

    aux=ftemp[m1++];
    fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
  }
};
// ****************************************************************************
void Linear::printtofile(FILE * fout,const double u[],
			 const Vertice vert[], const int map[])  // 10/02/2013
{
  double aux,x1,x2,x3;
  double eaux,aux1;
  double xa,xb;
  double ftemp[Q[0]];
  x2=0.0;
  x3=0.0;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;

  aux1=xb-xa;

  int m1=0;
  evalGQ(ftemp,u);
  for(int i=0;i<Q[0];i++){
    // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
    eaux=(xGQ[0][i]+1.0)/2.0;
    x1=eaux*aux1+xa;

    aux=ftemp[m1++];
    fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
  }
};

// ****************************************************************************
void Linear::printGQtofile(FILE * fout,const double ftemp[],
			     const double ftemp1[],
			     const Vertice vert[], const int map[])   // 10/02/2013
{
  double aux,x1,x2,x3;
  double eaux,aux1;
  double xa,xb;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  x2=0.0;
  x3=0.0;

  aux1=xb-xa;

  int m1=0;
  for(int i=0;i<Q[0];i++){
    // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
    eaux=(xGQ[0][i]+1.0)/2.0;
    x1=eaux*aux1+xa;
    aux=ftemp[m1];
    aux1=ftemp1[m1++];
    fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,aux1);
  }
};
// ****************************************************************************
// ****************************************************************************
void Linear::printwGQtofile(FILE * fout,
			   const Vertice vert[],
			    const int map[],
			    const double JV[])   // 10/02/2013
{
  double aux,x1,x2,x3;
  double eaux,aux1;
  double xa,xb;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  x2=0.0;
  x3=0.0;

  aux1=xb-xa;

  int m1=0;
  for(int i=0;i<Q[0];i++){
    // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
    eaux=(xGQ[0][i]+1.0)/2.0;
    x1=eaux*aux1+xa;
    aux=wGQ[0][i]*JV[m1++];
    fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
  }
};
// ****************************************************************************
void Linear::printtoarray(const double u[],
			    const Vertice vert[], const int map[],
			    double x[], double y[], double z[], double ftemp[])  // 10/02/2013
{
  //double x2,x3;
  double eaux,aux1;
  double xa,xb;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;

 // x2=0.0;
 // x3=0.0;
  aux1=xb-xa;

  int m1=0;
  evalGQ(ftemp,u);
  for(int i=0;i<Q[0];i++){
    // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
    eaux=(xGQ[0][i]+1.0)/2.0;
    /*x1*/x[m1]=eaux*aux1+xa;
    /*x2*/y[m1]=0.0;
    /*x3*/z[m1]=0.0;
    //aux=ftemp[m1];
    m1++;
    //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
  }
};
// ****************************************************************************
// Calcular os valores da funcao nos pontos de Gauss
// ****************************************************************************
void Linear::evalGQ(double f0[],double f1[],
			const double uh0[],const double uh1[])  // 10/02/2013
{
  double aux0,aux1,eta1,Fa;
  int a;
  int i,p;
  // ************************************************************************
  for(i=0;i<Q[0];i++){
    eta1=xGQ[0][i];
    aux0=0.0;
    aux1=0.0;

    // A
    a=0; p=0;
    Fa=Psia(P[0],p,eta1);
    aux0+=Fa*uh0[a];
    aux1+=Fa*uh1[a];
    // B
    a=1; p=P[0];
    Fa=Psia(P[0],p,eta1);
    aux0+=Fa*uh0[a];
    aux1+=Fa*uh1[a];
    // AB
    a=2;
    for(p=1;p<P[0];p++){
      Fa=Psia(P[0],p,eta1);
      aux0+=Fa*uh0[a];
      aux1+=Fa*uh1[a];
      a++;
    }
    f0[i]=aux0;
    f1[i]=aux1;
  }
};

// ****************************************************************************
// Calcular os valores da funcao nos nos de Gauss (overloaded)
// NF = numero de campos (default = 1)
// nvar= indice da variavel (default=0)
// ****************************************************************************
void Linear::evalGQ(double f0[],const double u0[],const int NF,
		      const int nvar)  // 10/02/2013
{
  double aux0,eta1,Fa;
  int a;
  int i,p;
  // ************************************************************************
  //cout << "EvalGQ:  u0 "<< u0[0] << " , " << u0[1] << "\n";
  for(i=0;i<Q[0];i++){
    eta1=xGQ[0][i];
    aux0=0.0;

    // A
    a=0; p=0;
    Fa=Psia(P[0],p,eta1);
    aux0+=Fa*u0[a*NF+nvar];
    //cout << "Fa(0) = " << Fa << "\n";
    // B
    a=1; p=P[0];
    Fa=Psia(P[0],p,eta1);
    aux0+=Fa*u0[a*NF+nvar];
    //cout << "Fa(P[0]) = " << Fa << "\n";
    // AB
    a=2;
    for(p=1;p<P[0];p++){
      Fa=Psia(P[0],p,eta1);
      aux0+=Fa*u0[a*NF+nvar];
      a++;
      //cout << "Fa(p) = " << Fa << "\n";
    }
    f0[i]=aux0;

    //cout << "i " << i << " eta1 "<< eta1 << " valor " << aux0 << "\n";

  }
};

// ****************************************************************************
// Evaluates the value of the field at the vertices of the element
// ****************************************************************************
void Linear::computeVertice(double f_vert[],const double u[],
		    const Vertice vert[], const int map[])  // 10/02/2013
{
  double aux, Fa,eta1;

  int m1;
  int i,p;
  for(i=0;i<2;i++){
    eta1=-1.0 + i*2.0;
    aux=0.0;
    for(m1=0;m1<nn;m1++){
      p= mode_[m1].p_val();
      Fa=Psia(P[0],p,eta1);
      aux+=Fa*u[m1];
    }
    //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
    f_vert[map[i]]=aux;
  }
};
// ****************************************************************************
// Evaluates the value of the field at points
// ****************************************************************************
void Linear::computeAtPoints(const int npoints, const double LocCoord[],const double u[],
			     const Vertice vert[], const int map[],double f[],double GloCoord[])
{
  double aux, Fa,eta1;

  int m1;
  int i,p;
  for(i=0;i<npoints;i++){
    eta1=LocCoord[i];
    aux=0.0;
    for(m1=0;m1<nn;m1++){
      p= mode_[m1].p_val();
      Fa=Psia(P[0],p,eta1);
      aux+=Fa*u[m1];
    }
    //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
    f[i]=aux;
  }
};
/*
void Linear::make_Phi(const int m,double Phi[])  // 10/02/2013
{
  double eta1;
  int i;
  int p;
  int n;
  p= mode_[m].p_val();
  n=0;
  for(i=0;i<Q[0];i++){
    eta1=xGQ[0][i];
    Phi[n++]=Psia(P[0],p,eta1);
  }
};
*/
void Linear::eval_Phi(const int n,double Phi[])  // 10/02/2013
{
  for(int i=0;i<NGQP;i++)Phi[i]=Phi_val[n][i];
};

void Linear::eval_GradPhi(const Vertice vert[], const int map[],const int m,double ** der)
{
};
// ****************************************************************************
void Linear::Jacobian(const Vertice vert[],const int map[],double * JV)  // 10/02/2013
{
  // Esta rotina eh valida so para linhas!
  double a1;
  double xa,xb;
  double Jacobian;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;

  a1=xb-xa;
  Jacobian=a1/2.0;
  if(Jacobian==0.0){
    printf("Jacobian=0\nVertices\n(%lf ), (%lf)\n",xa,xb);
  }
//  #ifdef PRINTF_ON
//    printf("Linear::Jacobian=%lf\n",Jacobian);
//  #endif
  for(int i=0;i<Q[0];i++)
    JV[i]=Jacobian;
};
/*
void Linear::Processar_geometria(int nel,
				 const Vertice * vert,
				 const int numv,
				 const int * VN,
				 int map[],
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
  int n0,n1;//ng0;
  int flag0,flag1;
  if(numv != 2){
    printf("Incompatibilidade: Dado do elemento nao eh de Segmento Linha\n");
    exit(0);
  }
  // passa o numero dos vertices
  n0=VN[0];
  flag0=0; // reset flag; flag = 0 equivale a aresta nao numerada
  for (int i=0; i<NL && flag0==0; i++){
    if(border[i].Na==n0){
      flag0=1;
      border[i].elemento[1]=nel;
      border[i].num_local[1]=0;
      border[i].sinal[1]=1;
      border[i].tipo=2;// aresta interior
    }
  }
  if(flag0==0){// aresta nova
    EDGE temp_border;

  //  border[NL].Na=n0;
  //  border[NL].Nb=n0;
  //  border[NL].elemento[0]=nel;
  //  border[NL].num_local[0]=0;
  //  border[NL].sinal[0]=1;
  //  border[NL].tipo=0;// aresta de contorno (inicialmente todas sao no-flow)

    temp_border.Na=n0;
    temp_border.Nb=n0;
    temp_border.elemento[0]=nel;
    temp_border.num_local[0]=0;
    temp_border.sinal[0]=1;
    temp_border.tipo=0;// aresta de contorno (inicialmente todas sao no-flow)
    border.push_back(temp_border);
    NL++;
    NG++;
  }

  n1=VN[1];
  flag1=0;
  for (int i=0; i<NL && flag1==0; i++){
    if(border[i].Na==n1){
      flag1=1;
      border[i].elemento[1]=nel;
      border[i].num_local[1]=1;
      border[i].sinal[1]=1;
      border[i].tipo=2;// aresta interior
    }
  }
  if(flag1==0){// aresta nova
    EDGE temp_border;

 //   border[NL].Na=n1;
//    border[NL].Nb=n1;
//    border[NL].elemento[0]=nel;
//    border[NL].num_local[0]=1;
//    border[NL].sinal[0]=1;
//    border[NL].tipo=0;// aresta de contorno (inicialmente todas sao no-flow)

    temp_border.Na=n1;
    temp_border.Nb=n1;
    temp_border.elemento[0]=nel;
    temp_border.num_local[0]=1;
    temp_border.sinal[0]=1;
    temp_border.tipo=0;// aresta de contorno (inicialmente todas sao no-flow)
    border.push_back(temp_border);
    NL++;
    NG++;
  }
  if(flag0*flag1==0) {
    // calculo da normal

    border[NL].comprimento=1.0;
    border[NL].normal[0] = 1.0;
    border[NL].normal[1] = 0.0;
  //  ng0=NG;
    // ajusta os modos dos vertices para fazer o gbnmap
    Ng[nel]=NG;
      sinal[0]=1;//sgn0;
    //make_gbnmap(n0,n1,ng0,map,sgn);
    NG+=(P[0]-1);// acrescenta numero de nos armazenados
  }
}
*/
// ****************************************************************************

// ****************************************************************************
/* OBSOLETA
void Linear::make_gbnmap(int n0,int n1,int ng0,int map[],int sgn[] )
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
  //printf("SINAIS de 0 1 2: %d %d %d\n", sgn[0],sgn[1],sgn[2]);

  imin=2;
  imax=p;
  // if(sign0==-1){
  temp=1;
  for(i=2; i<=p;i++){
    map[i]=ng0;
    sgn[i]=1;
    ng0++;
  }
};
*/
// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// ****************************************************************************

void Linear::Gradiente(FILE * fout, double * grad[],
			 const  double fvec[],
			 const Vertice vert[], const int map[])
{
  //cout << "Gradiente Linear "<< endl;
  double xa,xb,eta1,x1,e1p,aux1;
  double x2=0.0;
  double x3=0.0;

  //printf(" CALCULO DO GRADIENTE de um vetor. ELEMENTO LINEAR\n");

  Gradiente(grad,fvec,vert,map);

  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  aux1=(xb-xa);
  // Cheque do gradiente

  for(int i=0;i<Q[0];i++){
    eta1=xGQ[0][i];
    e1p=(1.0+eta1)/2.0;
    // coordenadas x1, x2, x3
    x1=xa+aux1*e1p;
    fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,grad[0][i], g(x1,x2,x3));
  }
};
// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// ****************************************************************************
// Necessita ser reescrito para considerar o caso do triangulo no espaco 3d

void Linear::Gradiente(double * grad[],
			 const  double fvec[],
			 const Vertice vert[], const int map[])
{
  //cout << "Gradiente Linear Principal "<< endl;
  double xa,xb;
  int l;
  double a1, J1D;
  double aux0;
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;

  // calculo dos coeficientes aij
  a1=(xb-xa);

  // calculo de J1D
  J1D=a1/2.0;

  // calculo das derivadas com relacao a eta1
  for(int p=0;p<Q[0];p++){
    aux0=0.0;
    for(l=0;l<Q[0];l++){
      aux0+=D[p][l][0]*fvec[l];
    }
    grad[0][p]=(aux0/J1D);
  }
};
// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[m][i] = component i of the gradient at the point m(pq)
// ****************************************************************************

void Linear::Gradiente(FILE * fout, double * grad[],
			 double (*func)(double, double, double),
			 const Vertice vert[], const int map[])
{
  //cout << "Gradiente 1 "<< endl;
  double xa,xb,eta1,x1;
  double x3=0.0;// triangulo no plano x,y
  double x2=0.0;
  double e1p;
  int i;
  double fvec[MAXQ];
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  double aux1= xb-xa;
  // calculo do vetor contendo a funcao nos pontos de integracao de Gauss

  for(i=0;i<Q[0];i++){
    eta1=xGQ[0][i];
    e1p=(1.0+eta1)/2.0;
    // coordenadas x1, x2, x3
    x1=(xa+aux1*e1p);
    //    x2=(ya+(yb-ya)*e1p);
    fvec[i]=func(x1,x2,x3);
  }
  // fprintf(fout," CALCULO DO GRADIENTE de uma funcao. ELEMENTO LINEAR\n");
  Gradiente(fout,grad,fvec,vert,map);
};

void Linear::print_nome(FILE * fout)
{
  fprintf(fout,"ELEMENTO LINEAR\n");
};
// ****************************************************************************
// Imposes Dirichlet type boundary conditions
// ****************************************************************************
void Linear::Dirichlet(const int no,
                       const Vertice vert[],
                       const int vert_map[],
                       const int nmap[],
                       const int sgn[],
                       int bflag[], // output
                       double X[],  // output
                       double (*f)(double,double,double))
{
  // **************************************************************************
  // flag = 0 : Dirichlet, valor conhecido, bflag=0
  //      = 1 : valor desconhecido, bflag=1
  // **************************************************************************
  int flag=0;
  int temp;
  double xa,ya,za;

  xa=vert[vert_map[no]].x;
  ya=vert[vert_map[no]].y;
  za=vert[vert_map[no]].z;

  temp=nmap[no];
  X[temp]=f(xa,ya,za);
  bflag[temp]=flag;
  cout << "Linear::Dirichlet: no " << no << " gbnmap " << temp << " valor " << X[temp]<< "\n";
};

void Linear::face_Jacobian(const int face_num,
                               const Vertice vert[],
                               const int vert_map[], // numero global dos vertices dos nos
                               const int sgn[],
                               double * J)
{};
void Linear::teste(int & v)
{
   v=100000;
};

// ****************************************************************************
// Calcula o traco e o coloca na ordem correta de acordo com o sinal da borda *
// ****************************************************************************
void Linear::trace(const int lado, const int qmax, const int sinal,
                   const double * valores,
                   double * saida,const int map[])
{
  if (lado ==1)saida[0] = valores[Q[0]];
  else saida[0] = valores[0];
};
// 10/02/2013

// ****************************************************************************
// Evaluates the value of the field at the Gauss Quadrature points
// ****************************************************************************
void Linear::computeFuncGQ(double f_[],
			   const Vertice vert[], const int map[],
			   double (*func)(double,double,double))  // 10/02/2013
{
  double x1,x2,x3;
  double eaux,aux1;
  double xa,xb;
  int i;//,count;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  x2=0.0;
  x3=0.0;
  aux1=xb-xa;
  for(i=0;i<Q[0];i++){
    eaux=(xGQ[0][i]+1.0)/2.0;
    // coordenadas (conforme Karniadakis pagina 108)
    x1=eaux*aux1+xa;
    f_[i]=func(x1,x2,x3);
  }
};
// *****************************************************************************
// Calcula os tracos de Phi, GradPhi e Jb nos pontos de Gauss sobre as arestas
// ordenando-os de forma que haja coincidencia de pontos dos elementos vizinhos
// *****************************************************************************
// Precisa verificar
void Linear::elem_traces(const Vertice vert[],const int map[],const int sinal[],
                         double *** TP,double **** TGP,double * Jb)
{
  // Precisa ser corrigido !!!!!!!!!!!!!!!!!!!!!!!!
  double xa,xb,eta1;
  int h,l,m;
    l=0;
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;

  double d1,aux0,aux1,der1;
 // int s0=nn*ndim;
  // int s1=   ndim;
  aux0=(xb-xa) / 2.0;
  for(h=0;h<nborder;h++){ // loop sobre as bordas
    switch(h) { // switch

    case 0:
      eta1=-1.0;
      break;

    case 1:
      eta1=1.0;
      break;

    }// switch


    for(m=0;m<nn;m++) { // loop sobre os modos
      int p= mode_[m].p_val();
    //  int i_m=h*nn+m;
      aux1=Psia(P[0],p,eta1);
      d1=DPsia(P[0],p,eta1);
      /*TP[i_m + l]*/TP[h][m][0]=aux1;
      der1 = d1/aux0;
      /*TGP[h*s0+m*s1+l]*/ TGP[h][nn][0][l] = der1;

      Jb[h+l]=aux0;
    } // loop sobre os modos
  } // loop sobre as bordas
}
// *******************************************
// Precisa verificar
void Linear::trace_Jb(const Vertice vert[],const int map[],const int sinal[],
                      double * Jb)
{
    // Precisa ser corrigido !!!!!!!!!!!!!!!!!!!!!!!!
    double xa,xb,eta1;
    int h;
    // coordenadas dos nos
    xa=vert[map[0]].x;
    xb=vert[map[1]].x;

    double d1,aux0,aux1,der1;
    // int s0=nn*ndim;
    // int s1=   ndim;
    aux0=(xb-xa) / 2.0;
    Jb[h]=aux0;

}
const int Linear::aresta_lvert(const int & i, const int & j) const {return aresta[i][j];};
const int Linear::face_lvert(const int & i, const int & j) const {return 0;};
const int Linear::show_nvf(const int &i) const {return 0;};
const int Linear::show_face_tipo(const int &i) const {return 0;};
const int Linear::show_fd0(const int &i) const {return 0;};
const int Linear::show_fd1(const int &i) const {return 0;};
const int Linear::show_fv2(const int &i) const {return 0;};

const int Linear::show_ind_mode(const int & i, const int & j, const int & k) const {return ind_mode_[i];};
void Linear::superficie_externa(const Vertice vert[],const int Vert_map[],
                                const int & num_local,
                                double & area,double normal[3])
{
  int v0,v1;
  int sinal_local =  sinal_normal[num_local];
  v0=Vert_map[0];
  v1=Vert_map[1];
  double lx,ly,lz;
  lx=vert[v1].x - vert[v0].x;
  ly=vert[v1].y - vert[v0].y;
  lz=vert[v1].z - vert[v0].z;

  normal[0] = sinal_local * lx;
  normal[1] = sinal_local * ly;
  normal[2] = sinal_local * lz;

  area = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

  normal[0]/=area;
  normal[1]/=area;
  normal[2]/=area;

};
//const int * Linear::show_face(const int &i){};
void Linear::face_GQCoord(const Vertice vert[],const int map[],
                          const int a0,const int qmax,
                          double x[],double y[],double z[])
{
};
