#include "spectral.h"
// ****************************************************************************
// Class Quadrilateral
// ****************************************************************************
// ****************************************************************************
Quadrilateral::Quadrilateral(int p0,int q0)
{
  //printf("Quadrilateral::Quadrilateral\n");
  tipo=3;
  ndim=2;
  nv=4;
  ne=4;
  nf=0;
  nborder=4;
  vtk_type=9;
  emapi = new int [ne+1];
  emapv = new int [ne*(p0+1)];
  bmapi = emapi;
  bmapv = emapv;
  set(p0, q0);
};
// ****************************************************************************
Quadrilateral::~Quadrilateral()
{
	//printf("Quadrilateral::~Quadrilateral\n");
  //printf("Destruindo Quadrilateral\n");
  delete [] emapi; emapi=nullptr;
  delete [] emapv; emapv=nullptr;
  bmapi=nullptr;
  bmapv=nullptr;
  delete [] D_Phi_val; D_Phi_val=nullptr;
  
  //libera memoria de ind_mode_
  for (int i=0;i<=P[0];++i){
    delete [] ind_mode_[i];  ind_mode_[i] = nullptr;
  }
  delete [] ind_mode_;  ind_mode_ = nullptr;
};
// ****************************************************************************

// Inicializa os dados do elemento
void Quadrilateral::set(int p0, int q0)
{
	//printf("Quadrilateral::set\n");
  int i,j,k;
  // Mapeamento das arestas
  // Indices que marcam o inicio daa arestas
  emapi[0]=0;
  emapi[1]=p0+1;
  emapi[2]=2*(p0+1);
  emapi[3]=3*(p0+1);
  emapi[4]=4*(p0+1);

  P[0]=p0; Q[0]=q0; gqt[0]=3;// Gauss-Lobatto-Jacobi
  P[1]=p0; Q[1]=q0; gqt[1]=3;// Gauss-Lobatto-Jacobi
  P[2]=1;  Q[2]=1;  gqt[2]=0;// Direcao nao usada
  NGQP=Q[0]*Q[1]*Q[2];
  qborder = q0;
  
  //aloca memoria para ind_mode_
  ind_mode_ = new int * [P[0]+1];
  for (int i=0;i<=P[0];++i){
    ind_mode_[i] = new int [P[1]+1];
  }

#ifdef PRINTF_ON
  printf("Em Quadrilateral::set(%d,%d)\n",p0,q0);
  printf("Em Quadrilateral::set apos Sair de Quadrilateral::setStdel\n");
  printf("P: %d %d\nQ: %d %d\ngqt : %d %d\n",P[0],P[1],Q[0],Q[1],gqt[0],gqt[1]);
#endif
  //    D ----- C
  //    |       |
  //    |       |
  //    A ----- B
  //
  mode_[0].set_mode(0,0,0);       ind_mode_[0][0] = 0;  // A
  mode_[1].set_mode(P[0],0,0);    ind_mode_[P[0]][0] = 1; // B
  mode_[2].set_mode(P[0],P[1],0); ind_mode_[P[0]][P[1]] = 2;// C
  mode_[3].set_mode(0,P[1],0);    ind_mode_[0][P[1]] = 3;// D
  int a=4;

  // AB
  k=emapi[0]; // aresta 0
  emapv[k]=0;// Vertice A
  k++;
  for(i=1;i<P[0];i++){
    mode_[a].set_mode(i,0,0);
    ind_mode_[i][0] = a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=1;// Vertice B

  // BC
  k=emapi[1];// aresta 1
  emapv[k]=1;// Vertice B
  k++;
  for(j=1;j<P[1];j++){
    mode_[a].set_mode(P[0],j,0);
    ind_mode_[P[0]][j] = a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=2;// Vertice C

  // DC
  k=emapi[2]; // aresta 2
  emapv[k]=3;// Vertice D
  k++;
  for(i=1;i<P[0];i++){
    mode_[a].set_mode(i,P[1],0);
    ind_mode_[i][P[1]] = a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=2;// Vertice C

  // AD
  k=emapi[3]; // aresta 3
  emapv[k]=0;// Vertice A
  k++;
  for(j=1;j<P[1];j++){
    mode_[a].set_mode(0,j,0);
    ind_mode_[0][j] = a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=3;// Vertice D

  nb=a;
  
  //Verificar
  // Interior modes : j runs fastest
  for(i=1; i<P[0]; i++){
    for(j=1; j<P[1];j++){
      mode_[a].set_mode(i,j,0);
      ind_mode_[i][j] = a;
      a++;
    }
  }
  
  nn=a;
  // if(nn==nmode)printf("OK! Calculo do nmode correto: %d, nn=%d\n",nmode,nn);
  // Calcula os parametros de Gauss
  gauss_parameters_default(); //alpha=beta=0.0
 
  if(nb>MAXNB ||nn>MAXMODES  ||nn-nb>MAXNI){
    printf("Ajuste parametros no spectral.h\n");
    printf("Numero de nos de fronteira do elemento (nb = %d)  (MAXNB = %d)\n",nb,MAXNB);
    printf("Numero de nos do elemento nn(= %d)                (MAXMODES  = %d)\n",nn,MAXMODES);
    printf("Numero de nos internos do elemento (ni = %d)      (MAXNI = %d)\n",nn-nb,MAXNI);
    exit(0);
  }
#ifdef PRINTF_ON
  printf("Em Quadrilateral::set  passou gauss_parameters\n");
  printf("nb=%d emapv tem %d valores\n",nb,k);
  printf("Saindo de Quadrilateral::set: nv = %d  nn = %d (< MAXMODES = %d) nb = %d (< MAXNB = %d)\n",nv,nn,MAXMODES,nb,MAXNB);
  printf("Quadrilateral::set nn = %d nb= %d\n",nn,nb);
#endif
 
  // Construcao da matriz Phi_val[nn][Q[0]*Q[1]]
  double eta1,eta2;
  double Fa[Q[0]],Fb[Q[1]];
  for(int m=0;m<nn;m++){
    
    int p= mode_[m].p_val();
    int q= mode_[m].q_val();
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      Fb[j]=Psia(P[1],q,eta2);
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
// 06/02/2008
// ****************************************************************************
void Quadrilateral::print(FILE * fout)
{
  printStdel();
	printf("Quadrilateral::print\n");
  fprintf(fout,"Quadrilateral::print\n");
  for(int i=0;i<nn_val(); i++){
    fprintf(fout,"no= %4d ",i);
    mode_[i].print(fout);
  }
  fprintf(fout,"\n");
  
};

// ****************************************************************************
double Quadrilateral::mass(int m1,int m2,const double JV[])
{
	//printf("Quadrilateral::mass(m1,m2,JV)\n");
  int p,q,r,s;
  //int flag1, flag2;
  int j, i;
  double aux, aaux, Fa,Fb,Ga,Gb;
   //printf("Entrou Quadrilateral::mass(%d, %d) JV\n",m1,m2);
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
      aaux+=wGQ[0][i]*Fa*Ga*JV[i+Q[0]*j];
    }
    Fb=Psia(P[1], q, xGQ[1][j]);
    Gb=Psia(P[1], s, xGQ[1][j]);
    aux+=wGQ[1][j]*Fb*Gb*aaux;
  }
  //printf("Saiu Quadrilateral::mass(%d, %d)\n",m1,m2);
  return aux;
};
//06/03/2008
// ****************************************************************************
// Creates the local matrices using static condensation
// It uses NEWMAT matrices and operators
// ****************************************************************************
void Quadrilateral::make_local_matrices()
{
	printf("Quadrilateral::make_local_matrices\n");
  int ni; 
  int i,ii,j,jj;
  ni=nn-nb;
#ifdef PRINTF_ON
  printf("Making local matrices 1 Quadrilateral: nb = %d  ni= %d\n", nb,ni);
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
	Fb=Psia(P[1],q, xGQ[1][j]);
	Gb=Psia(P[1],s, xGQ[1][j]);
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
  // newmat
#ifdef _NEWMAT
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
  // Matriz Mi
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
  printf("Saindo Quadrilateral::make_local_matrices: nb = %d  Ni= %d\n", nb,ni);
#endif
 
};

// ****************************************************************************

void Quadrilateral::make_mass_matrices(int NFields)
{
	printf("Quadrilateral::make_mass_matrices\n");
  /*
		// ***************************
    make_local_matrices();
    // ***************************
    if(NFields > 1)duplicar_mass_matrix(NFields);
    #ifdef PRINTF_ON
    printf("Saindo Quadrialteral::make_mass_matrices");
    printf("(%2d): nb = %d  ni= %d\n",NFields,nb,nn-nb);
    #endif
  */
};

//06/03/2008
// ****************************************************************************
// Estabelece os parametros de Gauss para o quadrilatero
// em coordenadas eta1 e eta2
// (Gauss-Jacobi quadratura) Ver Apendice B de Karniadakis & Sherwin
// ****************************************************************************
void Quadrilateral::gauss_parameters_default()
{
	//printf("Quadrilateral::gauss_parameters\n");
  double x[MAXQ];
  double wtemp[MAXQ];
  double Dtemp[MAXQ][MAXQ];

  // Valores para a terceira dimensao: 1 ponto (eta3=0.0) com w=1.0
  xGQ[2][0]=-1.0;
  wGQ[2][0]=1.0;
  // Valores para as duas dimensoes
  for(int i=0;i<ndim;++i){
   
    if(gqt[i]==1)
      {//gqt[i]=1 Gauss-Jacobi
        Gauss_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
      }
    
    else if(gqt[i]==2){
      // *******************************\/ \/\/\/\/****************************
      Gauss_Radau_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
    }
    
    else{//gqt[i]=3 
      Gauss_Lobatto_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
    }
    
    for(int j=0; j<Q[i]; ++j){
      xGQ[i][j]=x[j];
      wGQ[i][j]=wtemp[j];
    }
    
    for(int k=0;k<Q[i];++k){
      for(int l=0;l<Q[i];++l){
        D[k][l][i]=Dtemp[k][l];
	//printf(" Quadrilateral D[%d][%d][%d] = %g\n", k,l,i,D[k][l][i]);
      }
    }
  }
};
//06/03/2008
// ****************************************************************************
//
// ****************************************************************************
void Quadrilateral::vector_of_integral_of_f_Phi_dv(double vec[],
					 double (*func)(double,double,double),
					 const Vertice vert[],const int map[],
					 const double JV[])
{
	//printf("Quadrilateral::vector_of_integral_of_f_Phi_dv\n");
  int n,p,q;
  int j, i;
  double aux,eta1,eta2,x1,x2,x3;
  double a1m,a1p,a2m,a2p;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int Pdim=P[0]+1;
  double fp_xi2[Pdim][Q[1]];
  
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
 
  // construir matriz temporaria fp_xi2
  for(p=0;p<Pdim;p++){
    // fazer loop sobre xi2
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      a2m=(1.0-eta2)/2.0;
      a2p=(1.0+eta2)/2.0;
      //fazer soma em xi1
      aux=0.0;
      for(i=0;i<Q[0];i++){
	// coordenadas
	eta1=xGQ[0][i];
	a1m=(1-eta1)/2.0;
	a1p=(1+eta1)/2.0;
	x1= (xa*a1m + xb*a1p)*a2m + (xc*a1p + xd*a1m)*a2p;
	x2= (ya*a1m + yb*a1p)*a2m + (yc*a1p + yd*a1m)*a2p;
	x3= (za*a1m + zb*a1p)*a2m + (zc*a1p + zd*a1m)*a2p;
	aux+=(wGQ[0][i]*Psia(P[0],p,eta1)*func(x1,x2,x3)*JV[i+Q[0]*j]);
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
      aux+=(wGQ[1][j]*Psia(P[1], q, eta2)*fp_xi2[p][j]);
    }
    vec[n]=aux;
  }
};
//06/03/2008
// ****************************************************************************
//
// ****************************************************************************
void Quadrilateral::vector_of_integral_of_f_Phi_dv(double vec[],
					 const double func[],
					// const Vertice vert[],const int map[],
					 const double JV[])
{
	//printf("Quadrilateral::vector_of_integral_of_f_Phi_dv(func[])\n");
  int n,p,q;
  int j, i;
  double aux,eta1,eta2;//x1,x2,x3;
  int Pdim=P[0]+1;
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
	aux+=(wGQ[0][i]*Psia(P[0],p,eta1)*func[count]*JV[count]);
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
      aux+=(wGQ[1][j]*Psia(P[1], q, eta2)*fp_xi2[p][j]);
    }
    vec[n]=aux;
  }
};
//06/03/2008
// ****************************************************************************
void Quadrilateral::printtofile(FILE * fout,const double u[],
				double (*func)(double,double,double), 
				const Vertice vert[], const int map[])
{
	printf("Quadrilateral::printtofile\n");
  double aux, eta1,eta2,x1,x2,x3;
  double a1m,a1p,a2m,a2p;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  double ftemp[Q[0]*Q[1]];
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
  int m1=0;
  evalGQ(ftemp,u);
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    a2m=(1.0-eta2)/2.0;
    a2p=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
      eta1=xGQ[0][i];
      a1m=(1-eta1)/2.0;
      a1p=(1+eta1)/2.0;
      x1= (xa*a1m + xb*a1p)*a2m + (xc*a1p + xd*a1m)*a2p;
      x2= (ya*a1m + yb*a1p)*a2m + (yc*a1p + yd*a1m)*a2p;
      x3= (za*a1m + zb*a1p)*a2m + (zc*a1p + zd*a1m)*a2p;
      aux=ftemp[m1++];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
    }
  }
};
// ****************************************************************************
void Quadrilateral::printtofile(FILE * fout,const double u[],
				const Vertice vert[], const int map[])
{
	printf("Quadrilateral::printtofile1\n");
  double aux, eta1,eta2,x1,x2,x3;
  double eaux,faux,gaux,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  double ftemp[Q[0]*Q[1]];
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
  int m1=0;
  evalGQ(ftemp,u);
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/2.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
      eta1=xGQ[0][i];
      faux=(1-eta1)/2.0;
      gaux=(1+eta1)/2.0;
      x1=eaux*(xa*faux+xb*gaux)+haux*(xd*faux+xc*gaux);
      x2=eaux*(ya*faux+yb*gaux)+haux*(yd*faux+yc*gaux);
      x3=eaux*(za*faux+zb*gaux)+haux*(zd*faux+zc*gaux);
      aux=ftemp[m1++];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
    }
  }
};
// ****************************************************************************
void Quadrilateral::printGQtofile(FILE * fout,const double ftemp[],
				  const double ftemp1[],
				  const Vertice vert[], const int map[])
{
	//printf("Quadrilateral::printGQtofile\n");
  double aux, aux1,eta1,eta2,x1,x2,x3;
  double eaux,faux,gaux,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
  int m1=0;
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/2.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
      eta1=xGQ[0][i];
      faux=(1-eta1)/2.0;
      gaux=(1+eta1)/2.0;
      x1=eaux*(xa*faux+xb*gaux)+haux*(xd*faux+xc*gaux);
      x2=eaux*(ya*faux+yb*gaux)+haux*(yd*faux+yc*gaux);
      x3=eaux*(za*faux+zb*gaux)+haux*(zd*faux+zc*gaux);
      aux=ftemp[m1];
      aux1=ftemp1[m1++];
			
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,aux1);
    }
  }
};
// ****************************************************************************
// ****************************************************************************
void Quadrilateral::printwGQtofile(FILE * fout,
				   const Vertice vert[],
				   const int map[],
				   const double JV[])
{
	printf("Quadrilateral::printwGQtofile\n");
  double aux, eta1,eta2,x1,x2,x3;
  double eaux,faux,gaux,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
  int m1=0;
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/2.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
      eta1=xGQ[0][i];
      faux=(1-eta1)/2.0;
      gaux=(1+eta1)/2.0;
      x1=eaux*(xa*faux+xb*gaux)+haux*(xd*faux+xc*gaux);
      x2=eaux*(ya*faux+yb*gaux)+haux*(yd*faux+yc*gaux);
      x3=eaux*(za*faux+zb*gaux)+haux*(zd*faux+zc*gaux);
      aux=wGQ[0][i]*wGQ[1][j]*JV[m1++];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
    }
  }
};
// ****************************************************************************
void Quadrilateral::printtoarray(const double u[], 
			    const Vertice vert[], const int map[],
			    double x[], double y[], double z[], double ftemp[])
{
	printf("Quadrilateral::printtoarray\n");
  double eta1,eta2;
  double eaux,faux,gaux,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
  
  int m1=0;
  evalGQ(ftemp,u);
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/2.0;
    haux=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 108)
      eta1=xGQ[0][i];
      faux=(1-eta1)/2.0;
      gaux=(1+eta1)/2.0;
      /*x1*/x[m1]=eaux*(xa*faux+xb*gaux)+haux*(xd*faux+xc*gaux);
      /*x2*/y[m1]=eaux*(ya*faux+yb*gaux)+haux*(yd*faux+yc*gaux);
      /*x3*/z[m1]=eaux*(za*faux+zb*gaux)+haux*(zd*faux+zc*gaux);
      m1++;
      //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
    }
  }
};

// ****************************************************************************
// Calcular os valores da funcao nos pontos de Gauss
// ****************************************************************************
void Quadrilateral::evalGQ(double f0[],double f1[],
			const double uh0[],const double uh1[])
{
	printf("Quadrilateral::evalGQ0\n");
  double aux0,aux1, Fa,Fb,eta1,eta2;
  int a;
  int i,j,p,q,k;
  int n=0;
  int Pdim=P[0]+1;	
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
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*uh0[a];
    ftemp1[p]+=Fb*uh1[a];
    // B
    a=1;p=P[0];q=0;
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*uh0[a];
    ftemp1[p]+=Fb*uh1[a];
    // C 
    a=2;p=P[0];q=P[1]; 
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*uh0[a];
    ftemp1[p]+=Fb*uh1[a];
    // D 
    a=3; p=0; q=P[1]; 
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*uh0[a];
    ftemp1[p]+=Fb*uh1[a];    
    // AB
    a=4;
    q=0;
    Fb=Psia(P[1],q,eta2);
    for(p=1;p<P[0];p++){
      ftemp0[p]+=Fb*uh0[a];
      ftemp1[p]+=Fb*uh1[a];
      a++;
    }
    // BC
    p=P[0];
    for(q=1;q<P[1];q++){
      Fb=Psia(P[1],q,eta2);
      ftemp0[p]+=Fb*uh0[a];
      ftemp1[p]+=Fb*uh1[a];
      a++;
    }
    // DC
    q=P[1];
    Fb=Psia(P[1],q,eta2);
    for(p=1;p<P[0];p++){
      ftemp0[p]+=Fb*uh0[a];
      ftemp1[p]+=Fb*uh1[a];
      a++;
    }
    // AD
    p=0;
    for(q=1;q<P[1];q++){
      Fb=Psia(P[1],q,eta2);
      ftemp0[p]+=Fb*uh0[a];
      ftemp1[p]+=Fb*uh1[a];
      a++;
    }
    // Interior modes : j runs fastest
    for(p=1;p<P[0];p++){
        for(q=1;q<P[1];q++){
            Fb=Psia(P[1],q,eta2);
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
};
//06/03/2008
// ************************************************************************
// Calcular os valores da funcao nos nos de Gauss (overloaded)
// NF = numero de campos (default = 1)
// nvar= numero da variavel (default=0)
// ************************************************************************
void Quadrilateral::evalGQ(double f0[],const double u0[],const int NF,
			   const int nvar)
{
	//printf("Quadrilateral::evalGQ1\n");
  double aux0, Fa,Fb,eta1,eta2;
  int a;
  int i,j,p,q,k;
  int n=0;
  int Pdim=P[0]+1;
  double ftemp0[Pdim];
  //printf("Quadrilateral::evalGQ    NF = %d    nvar= %d\n",NF,nvar);
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j]; 
    // ********************************************************************
    // Construcao dos vetores temporarios
    // ********************************************************************
    for(k=0;k<Pdim;k++){
      ftemp0[k]=0.0;
    }
    // A
    a=0; p=0; q=0;
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*u0[a*NF+nvar];
    // B
    a=1; p=P[0]; q=0;
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*u0[a*NF+nvar];
    // C
    a=2; p=P[0]; q=P[1];
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*u0[a*NF+nvar];
  // D
    a=3;p=0;q=P[1];
    Fb=Psia(P[1],q,eta2);
    ftemp0[p]+=Fb*u0[a*NF+nvar];   
    // AB
    a=4;
    q=0; 
    Fb=Psia(P[1],q,eta2);
    for(p=1;p<P[0];p++){ 
      ftemp0[p]+=Fb*u0[a*NF+nvar];
      a++;
    }
    // BC
    p=P[0];
    for(q=1;q<P[1];q++){
      Fb=Psia(P[1],q,eta2);
      ftemp0[p]+=Fb*u0[a*NF+nvar];
      a++;
    }
    // CD
    q=P[1];
    Fb=Psia(P[1],q,eta2);
    for(p=1;p<P[0];p++){
      ftemp0[p]+=Fb*u0[a*NF+nvar];
      a++;
    }
    // AD
    p=0;
    for(q=1;q<P[1];q++){
      Fb=Psia(P[1],q,eta2);
      ftemp0[p]+=Fb*u0[a*NF+nvar];
      a++;
    }
    // Interior modes : j runs fastest
    for(p=1;p<P[0];p++){
      for(q=1;q<P[1];q++){
	Fb=Psia(P[1],q,eta2);
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
//06/03/2008
// ****************************************************************************
// Evaluates the value of the field at the vertices of the element
// ****************************************************************************
void Quadrilateral::computeVertice(double f_vert[],const double u[],
		    const Vertice vert[], const int map[])
{
	printf("Quadrilateral::computeVertice\n");
  double aux, Fa,Fb,eta1,eta2;//x1,x2,x3;
 // double eaux,faux,gaux,haux;
//  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int m1;
  int i,p,q;
 
//  xa=vert[map[0]].x;
//  xb=vert[map[1]].x;
//  xc=vert[map[2]].x;
//  xd=vert[map[3]].x;
//  ya=vert[map[0]].y;
//  yb=vert[map[1]].y;
//  yc=vert[map[2]].y;
//  yd=vert[map[3]].y;
//  za=vert[map[0]].z;
//  zb=vert[map[1]].z;
//  zc=vert[map[2]].z;
//  zd=vert[map[3]].z;
  
  i=0;
  for(eta2=-1.0;eta2<=1.0;eta2+=2.0){
//    eaux=(1.0-eta2)/2.0;
//    haux=(1.0+eta2)/2.0;
    for(eta1=eta2;eta1<=1.0;eta1+=2.0){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
//      faux=(1-eta1)/2.0;
//      gaux=(1+eta1)/2.0;
//      x1=eaux*(xa*faux+xb*gaux)+haux*(xd*faux+xc*gaux);
//      x2=eaux*(ya*faux+yb*gaux)+haux*(yd*faux+yc*gaux);
//      x3=eaux*(za*faux+zb*gaux)+haux*(zd*faux+zc*gaux);
      aux=0.0;
      for(m1=0;m1<nn;m1++){
	p= mode_[m1].p_val();
	q= mode_[m1].q_val();
	Fa=Psia(P[0],p,eta1);
	Fb=Psia(P[1],q,eta2);
        aux+=Fa*Fb*u[m1];
      }
      //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
      f_vert[map[i++]]=aux;
    }
  }
};
//06/03/2008
// ****************************************************************************
// Evaluates the value of the field at points
// ****************************************************************************
void Quadrilateral::computeAtPoints(const int npoints,
																		const double LocCoord[],
																		const double u[],
																		const Vertice vert[],
																		const int map[],
																		double f[],
																		double GloCoord[])
{
	printf("Quadrilateral::computeAtPoints\n");
  double aux, Fa,Fb,eta1,eta2,x1,x2,x3;
  double eaux,faux,gaux,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int m1;
  int i,p,q;
 
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
  
  for(i=0;i<npoints;i++){
    eta1=LocCoord[2*i ];
    eta2=LocCoord[2*i+1];
    eaux=(1.0-eta2)/2.0;
    haux=(1.0+eta2)/2.0;
    // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
    faux=(1-eta1)/2.0;
    gaux=(1+eta1)/2.0;
    x1=eaux*(xa*faux+xb*gaux)+haux*(xd*faux+xc*gaux);
    x2=eaux*(ya*faux+yb*gaux)+haux*(yd*faux+yc*gaux);
    x3=eaux*(za*faux+zb*gaux)+haux*(zd*faux+zc*gaux);
    GloCoord[3*i  ]=x1;
    GloCoord[3*i+1]=x2;
    GloCoord[3*i+2]=x3;
    aux=0.0;
    for(m1=0;m1<nn;m1++){
      p= mode_[m1].p_val();
      q= mode_[m1].q_val();
      Fa=Psia(P[0],p,eta1);
      Fb=Psia(P[1],q,eta2);
      aux+=Fa*Fb*u[m1];
    }
    //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
    //printf("%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
    f[i]=aux;
  }
};
// ****************************************************************************
// Evaluates the value of the field at the Gauss Quadrature points
// ****************************************************************************
void Quadrilateral::computeFuncGQ(double f_[], 
				  const Vertice vert[], const int map[],
				  double (*func)(double,double,double))
{
	//printf("Quadrilateral::computeFuncGQ\n");
  double eta1,eta2,x1,x2,x3;
  double eaux,faux,gaux,haux;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int i,j,count;
  
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;
  
  count=0;
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    eaux=(1.0-eta2)/2.0;
    haux=(1.0+eta2)/2.0;
    for(i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
      faux=(1-eta1)/2.0;
      gaux=(1+eta1)/2.0;
      x1=eaux*(xa*faux+xb*gaux)+haux*(xd*faux+xc*gaux);
      x2=eaux*(ya*faux+yb*gaux)+haux*(yd*faux+yc*gaux);
      x3=eaux*(za*faux+zb*gaux)+haux*(zd*faux+zc*gaux);
      f_[count++]=func(x1,x2,x3);
    }
  }
};
// ***************************************************************************/
//06/03/2008
/*
void Quadrilateral::make_Phi(const int m,double Phi[])
{
	printf("Quadrilateral::make_Phi\n");
  //printf("Quadrilateral::make_Phi, Q[0] = %d  Q[1]= %d, P[0] = %d, P[1] = %d\n",Q[0],Q[1],P[0],P[1]);
  double eta1,eta2;
  int p,q;
  int n;
  double Fa[Q[0]],Fb[Q[1]];
  p= mode_[m].p_val();
  q= mode_[m].q_val();
  //printf("AQUI\n");
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    Fb[j]=Psia(P[1],q,eta2);
  }
  //printf("AQUI\n");
  for(int i=0;i<Q[0];i++){
    eta1=xGQ[0][i];
    Fa[i]=Psia(P[0],p,eta1);
  }

  n=0;
  for(int j=0;j<Q[1];j++){
    for(int i=0;i<Q[0];i++){
      Phi[n++]=Fa[i]*Fb[j];
    }
  }
  //printf("saindo de make_Phi\n");
};
 */
//06/03/2008
void Quadrilateral::eval_Phi(const int n,double Phi[])
{
	//printf("Quadrilateral::eval_Phi\n");
  int q2=Q[0]*Q[1];
  for(int i=0;i<q2;i++)Phi[i]=Phi_val[n][i];
};
// ******************************************************************************
void Quadrilateral::eval_GradPhi(const Vertice vert[], const int map[],const int n,double ** der)
{
	printf("Quadrilateral::eval_GradPhi\n");
  double xa,ya,xb,yb,xc,yc,xd,yd,eta1,eta2;
 // double x3=0.0;
  double x1,x2,x12,y1,y2,y12;
  int m;
  double a11,a12,a21,a22, J2D;
  double b[2][2];
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  
  x1 =(-xa+xb+xc-xd)/4.0;    
  x12=( xa-xb+xc-xd)/4.0;
  x2 =(-xa-xb+xc+xd)/4.0;
  
  y1 =(-ya+yb+yc-yd)/4.0;
  y12=( ya-yb+yc-yd)/4.0;
  y2 =(-ya-yb+yc+yd)/4.0;

  int p = mode_[n].p_val();
  int q = mode_[n].q_val();
  double der1,der2,aux1,aux2,d1,d2;
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    a11=x1+x12*eta2;
    a21=y1+y12*eta2;
    aux2=Psia(P[1],q,eta2);
    d2= DPsia(P[1],q,eta2);
    for(int i=0;i<Q[0];i++){
      eta1=xGQ[0][i];
      aux1=Psia(P[0],p,eta1);
      d1= DPsia(P[0],p,eta1);
      der1 = d1*aux2;
      der2 = aux1*d2;
      a12=x2+x12*eta1;
      a22=y2+y12*eta1;
      J2D=a11*a22-a12*a21;
      m=i+j*Q[0];
      // calculo da matriz b[i][j]
      b[0][0]=a22/J2D;
      b[0][1]=-a12/J2D;
      b[1][0]=-a21/J2D;
      b[1][1]=a11/J2D;
      
      der[0][m] = der1 * b[0][0] + der2 * b[1][0];
      der[1][m] = der1 * b[0][1] + der2 * b[1][1];
    }
  } // loop sobre os pontos de Gauss
};

// ****************************************************************************


// ****************************************************************************
void Quadrilateral::Jacobian(const Vertice vert[],const int map[],double JV[])
{
  double x[3][4];
  double eta[2];
  double f[2][2]; // funcoes interpolantes
  double s[2] = { -0.5, 0.5}; // derivadas
  double dr[2][3];
  double aux,aux0,aux1;
  int a[2][4] = {
    {0,1,1,0},
    {0,0,1,1}
  };
  
  for(int i=0;i<4;++i)
  {
    x[0][i]=vert[map[i]].x;
    x[1][i]=vert[map[i]].y;
    x[2][i]=vert[map[i]].z;
  }
  
  int n=0;
  for(int j=0;j<Q[1];++j){
    eta[1]=xGQ[1][j];
    f[0][1] = 0.5*(1-eta[1]); // funcao inferior
    f[1][1] = 0.5*(1+eta[1]); // funcao superior
    for(int i=0;i<Q[0];++i){
      eta[0]=xGQ[0][i];
      f[0][0] = 0.5*(1-eta[0]); // funcao esquerda
      f[1][0] = 0.5*(1+eta[0]); // funcao direita
      
      for(int m=0;m<3;++m) {
        aux0=0.0;
        aux1=0.0;
        for(int l=0;l<4;++l){
          aux0 += x[m][l] * s[a[0][l]] * f[a[1][l]][1];
          aux1 += x[m][l] * f[a[0][l]][0] * s[a[1][l]];
        }
        dr[0][m]=aux0;
        dr[1][m]=aux1;
      }
  // Produto vetorial de dr[0][] por dr[1][]
      double v0 = dr[0][1] * dr[1][2] - dr[1][1] * dr[0][2];
      double v1 = dr[0][2] * dr[1][0] - dr[1][2] * dr[0][0];
      double v2 = dr[0][0] * dr[1][1] - dr[1][0] * dr[0][1];
      aux=sqrt(v0*v0 + v1*v1 + v2*v2);
      JV[n++] = aux;
     // cout << "Jacobiano = " << aux << endl;
      if(aux<=0.0)printf("Hexahedral::Jacobian Erro: Jacobiano nao-positivo = %g\n",aux);
   
    }
  }
};
// 01/09/2014


// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// ****************************************************************************

void Quadrilateral::Gradiente(double * grad[],
                              const  double fvec[],
                              const Vertice vert[], const int map[])
{
	//printf("Quadrilateral::Gradiente 1\n");
  double xa,ya,xb,yb,xc,yc,xd,yd,eta1,eta2;
 // double x3=0.0;
 // double a3p, a3m,a2p,a2m,a1p,a1m;
  double x1,x2,x12,y1,y2,y12;
  int i,j,l,m;
  double df[MAXQ*MAXQ][3];
  double a11,a12,a21,a22, J2D;
  double b[2][2];
  double aux0,aux1;
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;

  x1 =(-xa+xb+xc-xd)/4.0;
  x12=( xa-xb+xc-xd)/4.0;
  x2 =(-xa-xb+xc+xd)/4.0;

  y1 =(-ya+yb+yc-yd)/4.0;
  y12=( ya-yb+yc-yd)/4.0;
  y2 =(-ya-yb+yc+yd)/4.0;

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
  }
 
  // calculo das derivadas com relacao a x1 e x2
  for(int q=0;q<Q[1];q++){
    eta2=xGQ[1][q];
    a11=x1+x12*eta2;
    a21=y1+y12*eta2;
    for(int p=0;p<Q[0];p++){
      eta1=xGQ[0][p];
      a12=x2+x12*eta1;
      a22=y2+y12*eta1;
      J2D=a11*a22-a12*a21;
      m=p+q*Q[0];
      // calculo da matriz b[i][j]
      b[0][0]=a22/J2D;
      b[0][1]=-a12/J2D;
      b[1][0]=-a21/J2D;
      b[1][1]=a11/J2D;
      
      for(i=0;i<2;i++){
        aux0=0;
        for(j=0;j<2;j++){
          aux0+=b[j][i]*df[m][j];
        }
        grad[i][m]=aux0;
      }
    }
  }
};

// ****************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// ****************************************************************************

void Quadrilateral::Gradiente(FILE * fout, double * grad[],
			      const  double fvec[], 
			      const Vertice vert[], const int map[])
{
	printf("Quadrilateral::Gradiente 2\n");
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,eta1,eta2;
  double x3=0.0;
  double a2p,a2m,a1p,a1m;
  double x1,x2;//x12,y1,y2,y12;
  int m;
  
  // printf("Calculo do Gradiente de um vetor. Elemento Quadrilateral\n");
  Gradiente(grad,fvec,vert,map);
  
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  za=vert[map[0]].z;
  zb=vert[map[1]].z;
  zc=vert[map[2]].z;
  zd=vert[map[3]].z;

  // Cheque do gradiente
  
  for(int j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    a2m=(1.0-eta2)/2.0;
    a2p=(1.0+eta2)/2.0;
    for(int i=0;i<Q[0];i++){
      // coordenadas (conforme Karniadakis & Sherwin, pagina 109)
      eta1=xGQ[0][i];
      a1m=(1-eta1)/2.0;
      a1p=(1+eta1)/2.0;
      x1= (xa*a1m + xb*a1p)*a2m + (xc*a1p + xd*a1m)*a2p;
      x2= (ya*a1m + yb*a1p)*a2m + (yc*a1p + yd*a1m)*a2p;
      x3= (za*a1m + zb*a1p)*a2m + (zc*a1p + zd*a1m)*a2p;
      
      m=i+j*Q[0];
      fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,grad[0][m],grad[1][m],g(x1,x2,x3));
    }
  }
};

void Quadrilateral::Gradiente(FILE * fout, double * grad[],
			      double (*func)(double, double, double), 
			      const Vertice vert[], const int map[])
{
	printf("Quadrilateral::Gradiente 3\n");
  double xa,ya,xb,yb,xc,yc,xd,yd,eta1,eta2,x1,x2;
  double x3=0.0;// triangulo no plano x,y
  double e2p,e2m,e1p,e1m;
  int i,j,m;
  double fvec[MAXQ*MAXQ];
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
//  za=vert[map[0]].z;
//  zb=vert[map[1]].z;
//  zc=vert[map[2]].z;
//  zd=vert[map[3]].z;
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
      x1=(xa*e1m+xb*e1p)*e2m+e2p*(xc*e1p+xd*e1m);
      x2=(ya*e1m+yb*e1p)*e2m+e2p*(yc*e1p+yd*e1m);    
      m=i+j*Q[0];
      fvec[m]=func(x1,x2,x3);
    }
  }
  //printf("Calculo do Gradiente de uma funcao. Elemento Quadrilateral\n");
  Gradiente(fout,grad,fvec,vert,map);
};

void Quadrilateral::print_nome(FILE * fout)
{
	printf("Quadrilateral::print_nome\n");
  fprintf(fout,"ELEMENTO QUADRILATERAL\n");
}
// ****************************************************************************
// 20/03/2008

// ****************************************************************************
void Quadrilateral::Dirichlet(const int aresta,
                              const Vertice vert[],
                              const int vert_map[],
                              const int  nmap[],
                              const int sgn[],
                              int bflag[],
                              double X[],
                              double (*f)(double,double,double))
{
	//printf("Quadrilateral::Dirichlet\n");
  //cout << "Entrou em Quadrilateral::Dirichlet\n";
  int flag=0;
  // **************************************************************************
  // flag = 0 : Dirichlet, valor conhecido, bflag=0
  //      = 1 : valor desconhecido, bflag=1
  // **************************************************************************
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
  else if(aresta==2)// aresta 2; nos 3 e 2
    {
      xa=vert[vert_map[3]].x;
      xb=vert[vert_map[2]].x;
      ya=vert[vert_map[3]].y;
      yb=vert[vert_map[2]].y;
      za=vert[vert_map[3]].z;
      zb=vert[vert_map[2]].z;
      nd=0;//direcao 0
    }
  else // aresta 3; nos 0 e 3
    {
      xa=vert[vert_map[0]].x;
      xb=vert[vert_map[3]].x;
      ya=vert[vert_map[0]].y;
      yb=vert[vert_map[3]].y;
      za=vert[vert_map[0]].z;
      zb=vert[vert_map[3]].z;
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
  // valores nos vertices
  a0=f(xa,ya,za);
  ap=f(xb,yb,zb);

  if(p>1){// p>1
    gaux1=xa-xb;
    gaux2=ya-yb;
    gaux3=za-zb;
    J=sqrt(gaux1*gaux1+gaux2*gaux2+gaux3*gaux3)/2.0;
    
    for(k=0;k<q;k++)
      {
        eta1= x[k];//xGQ[0][k];
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
          for(k=0;k<q;k++) aux+=psi[i][k]*psi[j][k]*wtemp[k];//wGQ[0][k];
          aux*=J;
	
//	Ti[count]=i-1;
//	Tj[count]=j-1;
//	Tx[count++]=aux;
          A.element(i-1,j-1)=aux;
	
          if(j!=i){
//	  Ti[count]=j-1;
//	  Tj[count]=i-1;
//	  Tx[count++]=aux;
            A.element(j-1,i-1)=aux;
          }
        }
      
        aux=0.0;
        for(k=0;k<q;k++)
          aux+=psi[i][k]*(func[k]-a0*psi[0][k]-ap*psi[p][k])*wtemp[k];//wGQ[0][k];
        B.element(i-1)=aux*J;
      }
    
#ifdef _NEWMAT
      Y = A.i() * B; B=Y;
//#else 
//    ResolveSistema(p1,count,Ti,Tj,Tx,B,Y);
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
      ii=i+3;
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
      ii=i+P[0]+2;
      //temp=gbnmap[ii]*NFields+varn;
      temp=nmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  else if(aresta==2){
    //temp=gbnmap[3]*NFields+varn;
    temp=nmap[3];
    X[temp]=sgn[3]*a0;
    bflag[temp]=flag;
    //temp=gbnmap[2]*NFields+varn;
    temp=nmap[2];
    X[temp]=sgn[2]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++) {
      ii=i+P[0]+P[1]+1;
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
    //temp=gbnmap[3]*NFields+varn;
    temp=nmap[3];
    X[temp]=sgn[3]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++) {
      ii=i+2*P[0]+P[1];
      //temp=gbnmap[ii]*NFields+varn;
      temp=nmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  //cout << "Terminou Quadrilateral::Dirichlet\n";
};
// 20/03/2008

void Quadrilateral::teste(int & v)
{
	printf("Quadrilateral::teste\n");
   v=100000;
}

//Verificar
/*
// ****************************************************************************
// Mass Matrix entry M(m1, m2)
// ****************************************************************************
// Evaluates the inner product of Phi_m1 by Phi_m2
// ****************************************************************************
double Quadrilateral::mass(int m1,int m2)
{
 //printf("Quadrilateral::mass\n");
  int p,q,r,s;
  int flag1, flag2;
  int j, i;
  double aux, aaux, Fa,Fb,Ga,Gb;
  // printf("\nEntrou Quadrilateral::mass(%d, %d)\n",m1,m2);
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
      aaux+=wGQ[0][i]*Fa*Ga;// *JV[i+Q[0]*j];// <==-multiplica pelo Jacobiano-< 
    }
    Fb=Psia(P[1], q, xGQ[1][j]);
    Gb=Psia(P[1], s, xGQ[1][j]);
    aux+=wGQ[1][j]*Fb*Gb*aaux;
  }
  //printf("Saiu Quadrilateral::mass(%d, %d)\n",m1,m2);
 return (aux*JV[0]);// <----------multiplica pelo Jacobiano----<
}; */
//06/03/2008
// *****************************************************************************
// Calcula os tracos de Phi, GradPhi e Jb nos pontos de Gauss sobre as arestas
// ordenando-os de forma que haja coincidencia de pontos dos elementos vizinhos
// *****************************************************************************
void Quadrilateral::elem_traces(const Vertice vert[],const int map[],const int sinal[],
				double *** TP,double **** TGP,
				double * Jb)
{
	//printf("Quadrilateral::elem_traces\n");
  double xa,ya,xb,yb,xc,yc,xd,yd,eta1[qborder],eta2[qborder];
 // double x3=0.0;
  double x1,x2,x12,y1,y2,y12;
  int h,i,l,m;
  double a11,a12,a21,a22, J2D;
  double b[2][2];
  // coordenadas dos nos
  xa=vert[map[0]].x;
  xb=vert[map[1]].x;
  xc=vert[map[2]].x;
  xd=vert[map[3]].x;
  ya=vert[map[0]].y;
  yb=vert[map[1]].y;
  yc=vert[map[2]].y;
  yd=vert[map[3]].y;
  
  x1 =(-xa+xb+xc-xd)/4.0;    
  x12=( xa-xb+xc-xd)/4.0;
  x2 =(-xa-xb+xc+xd)/4.0;
  
  y1 =(-ya+yb+yc-yd)/4.0;
  y12=( ya-yb+yc-yd)/4.0;
  y2 =(-ya-yb+yc+yd)/4.0;
  
  //pontos de Gauss em uma dimensao
  double x[qborder], wtemp[qborder], Dtemp[MAXQ][MAXQ];
  Gauss_Jacobi_parameters(qborder, 0.0, 0.0, x, wtemp, Dtemp);
  
  double d1,d2,aux0,aux1,aux2,der1,der2;
 // int s0=nn*ndim*qborder;
  //int s1=   ndim*qborder;
  
  for(h=0;h<nborder;h++){ // loop sobre as bordas
    switch(h) { // switch
    
    case 0:
      aux0=sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya)) / 2.0;
      for(i=0;i<qborder;i++){
				eta1[i]=x[i]*sinal[h];
				eta2[i]=-1.0;
      }
      break;

    case 1:
      aux0=sqrt( (xb-xc)*(xb-xc) + (yb-yc)*(yb-yc)) / 2.0;
      for(i=0;i<qborder;i++){
				eta1[i]=1.0;
				eta2[i]=x[i]*sinal[h];
      }
      break;
    
    case 2:
      aux0=sqrt( (xc-xd)*(xc-xd) + (yc-yd)*(yc-yd)) / 2.0;
      for(i=0;i<qborder;i++){
				eta1[i]=x[i]*sinal[h];
				eta2[i]=1.0;
      }
      break;

    case 3:
      aux0=sqrt( (xd-xa)*(xd-xa) + (yd-ya)*(yd-ya)) / 2.0;
      for(i=0;i<qborder;i++){
				eta1[i]=-1.0;
				eta2[i]=x[i]*sinal[h];
      }
      break;

    }// switch

    for(l=0;l<qborder;l++){ // loop sobre os pontos de Gauss
      Jb[h*qborder+l]=aux0;
    }

    for(m=0;m<nn;m++) { // loop sobre os modos
      int p= mode_[m].p_val();
      int q= mode_[m].q_val();

    //  int i_m=h*nn*qborder+m*qborder;

      for(l=0;l<qborder;l++){ // loop sobre os pontos de Gauss
				aux1=Psia(P[0],p,eta1[l]);
				d1= DPsia(P[0],p,eta1[l]);
				aux2=Psia(P[1],q,eta2[l]);
				d2= DPsia(P[1],q,eta2[l]);
				TP[h][m][l]=aux1*aux2;
				//TP[i_m + l]=aux1*aux2;
				der1 = d1*aux2;
				der2 = aux1*d2;
	
				a11=x1+x12*eta2[l];
				a21=y1+y12*eta2[l];
				a12=x2+x12*eta1[l];
				a22=y2+y12*eta1[l];
				J2D=a11*a22-a12*a21;
				// calculo da matriz b[i][j]
				b[0][0]=a22/J2D;
				b[0][1]=-a12/J2D;
				b[1][0]=-a21/J2D;
				b[1][1]=a11/J2D;
				TGP[h][m][0][l] = der1 * b[0][0] + der2 * b[1][0];
				TGP[h][m][1][l] = der1 * b[0][1] + der2 * b[1][1];
	//TGP[h*s0+m*s1     +l] = der1 * b[0][0] + der2 * b[1][0];
	//TGP[h*s0+m*s1+qborder+l] = der1 * b[0][1] + der2 * b[1][1];
      } // loop sobre os pontos de Gauss
      
    } // loop sobre os modos

  } // loop sobre as bordas
}
// ****************************************************************************
void Quadrilateral::trace_Jb(const Vertice vert[],const int map[],const int sinal[],
                             double * Jb)
{
    //printf("Quadrilateral::elem_traces\n");
    double xa,ya,xb,yb,xc,yc,xd,yd,eta1[qborder],eta2[qborder];
    // double x3=0.0;
    double x1,x2,x12,y1,y2,y12;
    int h,i,l,m;
    double a11,a12,a21,a22, J2D;
    double b[2][2];
    // coordenadas dos nos
    xa=vert[map[0]].x;
    xb=vert[map[1]].x;
    xc=vert[map[2]].x;
    xd=vert[map[3]].x;
    ya=vert[map[0]].y;
    yb=vert[map[1]].y;
    yc=vert[map[2]].y;
    yd=vert[map[3]].y;
    
    x1 =(-xa+xb+xc-xd)/4.0;
    x12=( xa-xb+xc-xd)/4.0;
    x2 =(-xa-xb+xc+xd)/4.0;
    
    y1 =(-ya+yb+yc-yd)/4.0;
    y12=( ya-yb+yc-yd)/4.0;
    y2 =(-ya-yb+yc+yd)/4.0;
    
    //pontos de Gauss em uma dimensao
    double x[qborder], wtemp[qborder], Dtemp[MAXQ][MAXQ];
    Gauss_Jacobi_parameters(qborder, 0.0, 0.0, x, wtemp, Dtemp);
    
    double d1,d2,aux0,aux1,aux2,der1,der2;
    // int s0=nn*ndim*qborder;
    //int s1=   ndim*qborder;
    
    for(h=0;h<nborder;h++){ // loop sobre as bordas
        switch(h) { // switch
                
            case 0:
                aux0=sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya)) / 2.0;
                for(i=0;i<qborder;i++){
                    eta1[i]=x[i]*sinal[h];
                    eta2[i]=-1.0;
                }
                break;
                
            case 1:
                aux0=sqrt( (xb-xc)*(xb-xc) + (yb-yc)*(yb-yc)) / 2.0;
                for(i=0;i<qborder;i++){
                    eta1[i]=1.0;
                    eta2[i]=x[i]*sinal[h];
                }
                break;
                
            case 2:
                aux0=sqrt( (xc-xd)*(xc-xd) + (yc-yd)*(yc-yd)) / 2.0;
                for(i=0;i<qborder;i++){
                    eta1[i]=x[i]*sinal[h];
                    eta2[i]=1.0;
                }
                break;
                
            case 3:
                aux0=sqrt( (xd-xa)*(xd-xa) + (yd-ya)*(yd-ya)) / 2.0;
                for(i=0;i<qborder;i++){
                    eta1[i]=-1.0;
                    eta2[i]=x[i]*sinal[h];
                }
                break;
                
        }// switch
        
        for(l=0;l<qborder;l++){ // loop sobre os pontos de Gauss
            Jb[h*qborder+l]=aux0;
        }
    }
}
// ****************************************************************************
// Calcula o traco e o coloca na ordem correta de acordo com o sinal da borda *
// ****************************************************************************
void Quadrilateral::trace(const int lado,const int qmax,const int sinal,
                          const double *valores,double *saida)
{
    //printf("Quadrilateral::trace\n");
    int nd,ind,inc;

    if(lado==0 || lado==2){
        nd=0;
        }
    else nd=1;
  
    int q=Q[nd];
    double temp[q];
  
    if(lado==0){
        ind=0;
        inc=1;
    }
    else if (lado==1){
        ind=q-1;
        inc=q;
    }
    else if(lado==2){
        ind=q*(q-1);
        inc=1;
    }
    else{
        ind=0;
        inc=q;
    }
  
    // *****************************************
    // Para considerar o sinal de percurso das
    // arestas usa-se a expressao seguinte
    if(gqt[nd]==3){
        
       // if(sinal == -1){
       //     ind = ind + (q-1) * inc;
       //     inc = sinal * inc;
        //}
        ind = ind + (1 - sinal)/2 * (q-1) * inc;
        inc = inc * sinal;
    }
    // ******************************************
    // Transfere os dados de entrada para a saida
    for(int i=0;i<q;++i){
        saida[i]=valores[ind];
        ind+=inc;
    }
    // **********************************************
    // Se precisar fazer interpolação
    if(gqt[nd] != 3){
    // Interpolar para conter os pontos extremos
    // Executa quando os pontos de Gauss do elemento nao contem os extremos
        double Old[q];
        for(int i = 0; i < q; i++){
            temp[i]=saida[i];
            Old[i]=xGQ[nd][i];
        }
        for(int k=0;k<q;k++){
            double prod=1.0;
            double y=Old[k];
            for(int j=0;j<q;j++){
                if(j!=k)prod*=(y-Old[j]);
            }
            temp[k]/=prod;
        }
        // ****************************************************************
        // Os pontos de Gauss onde estao os tracos sao do tipo
        // Gauss-Jacobi (nao inclui os pontos dos extremos),
        // portanto sao diferentes dos pontos de Gauss
        // do elemento. Logo toda aresta necessita calcular o traco
        // ****************************************************************
        double Jac[qmax],wtemp[qmax],Dtemp[MAXQ][MAXQ];
        Gauss_Jacobi_parameters(qmax,0.0,0.0,Jac,wtemp,Dtemp);
    
        for(int i=0;i<qmax;i++){
            double sum=0.0;
            double y=Jac[i] * sinal; // percorre os pontos de Gauss na ordem decrescente se o sinal == -1
            for(int k=0;k<q;k++){
                double prod=1.0;
                for(int j=0;j<q;j++){
                    if(j!=k)prod*=(y-Old[j]);
                }
                sum+=(temp[k]*prod);
            }
            saida[i]=sum;
        }
    }
};
// revisado em 25/10/2011
const int Quadrilateral::aresta_lvert(const int & i, const int & j) const {return aresta[i][j];};
const int Quadrilateral::face_lvert(const int & i, const int & j) const {return face[i][j];};
const int Quadrilateral::show_nvf(const int &i) const {return 0;};
const int Quadrilateral::show_face_tipo(const int &i) const {return 0;};
const int Quadrilateral::show_fd0(const int &i) const {return 0;};
const int Quadrilateral::show_fd1(const int &i) const {return 0;};
const int Quadrilateral::show_fv2(const int &i) const {return 0;};
const int Quadrilateral::show_ind_mode(const int & i, const int & j, const int & k) const {return ind_mode_[i][j];};
void Quadrilateral::superficie_externa(const int Vert_map[],const Vertice vert[],
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
