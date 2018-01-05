#include "spectral.h"
// ***************************************************************************
// Class Tetraedro
// ***************************************************************************
Tetrahedral::Tetrahedral(int P, int Q)
{
  tipo=4;
  ndim=3;
  nv=4;
  ne=6;
  nf=4;
  nborder=4; // = nf
  vtk_type=10;
	emapi = new int [ne+1];
  emapv = new int [ne*(P+1)];
  set(P,Q);
  // Stdel::set_Typeid(vtk_type);
};

Tetrahedral::~Tetrahedral()
{
  //printf("Tetrahedral::~Tetrahedral\n");
  //printf("Destruindo Tetrahedral\n");
  delete [] emapi; emapi=nullptr;
  delete [] emapv; emapv=nullptr;
  delete [] D_Phi_val; D_Phi_val=nullptr;
  
  //libera memoria de ind_mode_
  for (int i=0;i<=(P[0]+1);++i){
    for(int j=0;j<=(P[1]+1);++j){
      delete [] ind_mode_[i][j];  ind_mode_[i][j] = nullptr;
    }
    delete [] ind_mode_[i];  ind_mode_[i] = nullptr;
  }
  delete [] ind_mode_;  ind_mode_ = nullptr;

  //Stdel::~Stdel();
};

void Tetrahedral::set(int p, int q)
{
  // Mapeamento das arestas
  // Indices que marcam o inicio das arestas
  emapi[0]=0;
  emapi[1]=p+1;
  emapi[2]=2*(p+1);
  emapi[3]=3*(p+1);
  emapi[4]=4*(p+1);
  emapi[5]=5*(p+1);
  emapi[6]=6*(p+1);
  
  printf("Em Tetrahedral::set\n");
  P[0]=p; Q[0]=q; gqt[0]=3;//Gauss-Lobatto
  P[1]=p; Q[1]=q; gqt[1]=2;//Gauss-Radau
  P[2]=p; Q[2]=q; gqt[2]=2;//Gauss-Radau
  NGQP=q*q*q;
  qborder = q*q;
  
 //aloca memoria para ind_mode_
  ind_mode_ = new int ** [P[0]+2];
  for (int i=0;i<=(P[0]+1);++i){
    ind_mode_[i] = new int * [P[1]+2];
    for(int j=0;j<=(P[1]+1);++j){
      ind_mode_[i][j] = new int [P[2]+1];
    }
  }
#ifdef PRINTF_ON
  printf("Em Tetrahedral::set apos Sair de Tetrahedral::setStdel\n");
  printf("P: %d %d %d \nQ: %d %d %d\ngqt : %d %d %d\n",P[0],P[1],P[2],Q[0],Q[1],Q[2],gqt[0],gqt[1],gqt[2]);
#endif  
  mode_[0].set_mode(0,0,0); ind_mode_[0][0][0] = 0;// A
  mode_[1].set_mode(P[0],0,0); ind_mode_[P[0]][0][0]=1;// B
  mode_[2].set_mode(P[0]+1,P[1],0);ind_mode_[P[0]+1][P[1]][0]=2;// C ponto colapsado
  mode_[3].set_mode(P[0]+1,P[1]+1,P[2]);ind_mode_[P[0]+1][P[1]+1][P[2]]=3;// D   ponto colapsado
  int a=4;
  
  // ****************************************************************************************
  // Sequencia de arestas especificado em Tetrahedral.h (aresta[][], ad0[], av1[], av2[] )
  // ****************************************************************************************
  for(int i=0;i<5;++i) { // loop sobre as arestas exceto a ultima, que eh colapsada
    int k=emapi[i];
    emapv[k]=aresta[i][0];
    k++;
    int dir[3] = {ad0[i], (ad0[i]+1)%3, (ad0[i]+2)%3};
    int ind[3];
    ind[dir[1]] = av1[i] == 0 ? 0 : P[dir[1]];
    ind[dir[2]] = av2[i] == 0 ? 0 : P[dir[2]];
    for(int j=1;j<P[dir[0]];++j) {
      ind[dir[0]] = j;
      mode_[a].set_mode(ind[0],ind[1],ind[2]);
      ind_mode_[ind[0]][ind[1]][ind[2]]=a;
      emapv[k]=a;
      k++;
      a++;
    }
    emapv[k]=aresta[i][1];
  }
  // CD; aresta colapsada; necessita especificacao especial
  a_CD=a;
  int k=emapi[5]; // sexta aresta
  emapv[k]=2; // Vertice C
  k++;
  for(int l=1;l<P[2];l++){
    mode_[a].set_mode(P[0]+1,P[1],l);// <=======
    ind_mode_[P[0]+1][P[1]][l]=a; // <=======
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=3; // Vertice D
  
  // ***************************************************************************************
  /* // Sequencia fixa de arestas
   // **************************************************************************************
 // cout << "inicio das arestas a= "<< a<< endl;
  // AB
  int k=0; // primeira aresta
  emapv[k]=0; // ponto A
  k++;
  for(int i=1;i<P[0];i++){
    mode_[a].set_mode(i,0,0); ind_mode_[i][0][0]=a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=1; // Vertice B
  
  // AC
  k=emapi[2]; // segunda aresta
  emapv[k]=0; // Vertice A
  k++;
  for(int j=1;j<P[1];j++){
    mode_[a].set_mode(0,j,0); ind_mode_[0][j][0]=a;
    emapv[k]=a;
    k++;
    a++;
  }
   emapv[k]=2; // Vertice C
   //cout << "a3 = "<< a << endl;
  
  // AD
  k=emapi[3]; // terceira aresta
  emapv[k]=0; // Vertice A
  k++;
  for(int l=1;l<P[2];l++){
    mode_[a].set_mode(0,0,l); ind_mode_[0][0][l]=a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=3; // Vertice D
   //cout << "a4 = "<< a << endl;
  
  //cout << "a1 = "<< a << endl;
  // BC
  k=emapi[1]; // quarta  aresta
  emapv[k]=1; // Vertice B
  k++;
  for(int j=1;j<P[1];j++){
    mode_[a].set_mode(P[0],j,0);ind_mode_[P[0]][j][0]=a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=2; // Vertice C
  //cout << "a2 = "<< a << endl;
  
  // BD
  k=emapi[4]; // quinta aresta
  emapv[k]=1; // Vertice B
  k++;
  for(int l=1;l<P[2];l++){
    mode_[a].set_mode(P[0],0,l); ind_mode_[P[0]][0][l]=a;
    emapv[k]=a;
    k++;
    a++;
  }
   emapv[k]=3; // Vertice D
   //cout << "a5 = "<< a << endl;
  
  // CD; aresta colapsada
  a_CD=a;
  k=emapi[5]; // sexta aresta
  emapv[k]=2; // Vertice C
  k++;
  for(int l=1;l<P[2];l++){
    mode_[a].set_mode(P[0]+1,P[1],l); ind_mode_[P[0]+1][P[1]][l]=a;
    emapv[k]=a;
    k++;
    a++;
  }
  emapv[k]=3; // Vertice D
  */
  
  //cout << "terminou arestas a6 = "<< a << endl;
  
  // Surface modes : r runs fastest
  //ABC Face 0
   for(int i=1; i<P[0]; i++){
    for(int j=1; j<P[1]-i;j++){
      mode_[a].set_mode(i,j,0); ind_mode_[i][j][0]=a;
      a++;
    }
  }
  //cout << "a7 = "<< a << endl;
  // ABD Face 1
  for(int i=1; i<P[0]; i++){
    for(int k=1; k<P[2]-i;k++){
      mode_[a].set_mode(i,0,k);ind_mode_[i][0][k]=a;
      a++;
    }
  }
  //cout << "a9 = "<< a << endl;
  // ACD Face 3 -> 2
  for(int j=1; j<P[1]; j++){
    for(int k=1; k<P[2]-j;k++){
      mode_[a].set_mode(0,j,k); ind_mode_[0][j][k]=a;
      a++;
    }
  }
  //cout << "a8 = "<< a << endl;
  //BCD Face 2 -> 3
  for(int j=1; j<P[1]; j++){
    for(int k=1; k<P[2]-j;k++){
      mode_[a].set_mode(P[0],j,k); ind_mode_[P[0]][j][k]=a;
      a++;
    }
  }
 
 nb=a;
  //cout << "a10 = "<< a << endl;
  // Interior modes: r runs fastest
  for(int i=1; i<P[0]; i++){
    for(int j=1; j<P[1]-i;j++){
      for(int k=1;k<P[2]-i-j;k++){
	mode_[a].set_mode(i,j,k);
        ind_mode_[i][j][k]=a;
	a++;
      }
    }
  }
  nn=a;
  //cout << "a11 = "<< a << endl;
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
  printf("Saindo de Tetrahedral::set: nv = %d  nn = %d (< MAXMODES = %d) nb = %d (< MAXNB = %d)\n",nv,nn,MAXMODES,nb,MAXNB);
#endif
};

// ************************************************************************
// ************************************************************************
// Sets Gauss parameters for the tetrahedral element
// Ja em coordenadas eta1,eta2 e eta3, 
// incluindo os termo (1-eta2)/2 e ((1-eta3)/2)^2 
// (Gauss-Jacobi quadratura) Ver Apendice B de Karniadakis & Sherwin
// ************************************************************************
void Tetrahedral::gauss_parameters_default()
{
   double x[MAXQ];
   double wtemp[MAXQ];
   double Dtemp[MAXQ][MAXQ];

   // Primeira coordenada
   Gauss_Lobatto_Jacobi_parameters(Q[0],0.0,0.0,x,wtemp,Dtemp);
   for(int j=0; j<Q[0]; ++j){
     xGQ[0][j]=x[j];
     wGQ[0][j]=wtemp[j];
   }
   for(int k=0;k<Q[0];++k)
     for(int l=0;l<Q[0];++l){
       D[k][l][0]=Dtemp[k][l];
     }   
   // Segunda coordenada
   // Ja em coordenadas eta2, incluindo o termo (1-eta2)/2
   // alpha=1.0, beta=0.0            \/ \/\/\/\/***************************
   Gauss_Radau_Jacobi_parameters(Q[1],1.0,0.0,x,wtemp,Dtemp);
   // ********************************^^^^^^^^^****************************
   for(int j=0; j<Q[1]; ++j){
     xGQ[1][j]=x[j];
     // *****\/\/\/\/\/\/**************************************************
     wGQ[1][j]=wtemp[j]/2.0;
     // *****^^^^^^^^^^^^**************************************************
   }
   for(int k=0;k<Q[1];++k)
     for(int l=0;l<Q[1];++l){
       D[k][l][1]=Dtemp[k][l];
     }
   // Terceira coordenada
   // Ja em coordenada eta3, incluindo o termo ((1-eta3)/2)^2
   // alpha=2.0, beta=0.0            *\/ \/\/\/\/**************************
   Gauss_Radau_Jacobi_parameters(Q[2],2.0,0.0,x,wtemp,Dtemp);
   // ********************************^^^^^^^^^****************************
   for(int j=0; j<Q[2]; ++j){
     xGQ[2][j]=x[j];
     // *****\/\/\/\/\/\/**************************************************
     wGQ[2][j]=wtemp[j]/4.0;
     // *****^^^^^^^^^^^^**************************************************
   }
   for(int k=0;k<Q[2];++k)
     for(int l=0;l<Q[2];++l){
       D[k][l][2]=Dtemp[k][l];
     } 
};
// ************************************************************************
// ************************************************************************
void Tetrahedral::vector_of_integral_of_f_Phi_dv(double vec[],
				       double (*func)(double,double,double), 
				       const Vertice vert[], const int map[],
				       const double JV[])
{
  int n,p,q,r;
  int i,j,k;
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  double e1p,e1m,e2p,e2m,e3p,e3m;
  double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  double xa,ya,xb,yb,xc,yc,xd,yd,za,zb,zc,zd;
  int Pdim=P[0]+2;
  double fp[Pdim][Q[1]][Q[2]], fpq[Pdim][Pdim][Q[2]];
  
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
  x1=xd-xa;
  x2=xb-xa;
  x3=xc-xa;
  y1=yd-ya;
  y2=yb-ya;
  y3=yc-ya;
  z1=zd-za;
  z2=zb-za;
  z3=zc-za;
 
  // primeira matriz fp(xi2,xi3)
  for(p=0;p<Pdim;p++){
    // fazer loop sobre xi2 e xi3 
    for(k=0;k<Q[2];k++){
      eta3=xGQ[2][k];
      e3p=(1.0+eta3)/2.0;
      e3m=(1.0-eta3)/2.0;
      for(j=0;j<Q[1];j++){
	eta2=xGQ[1][j];
	e2p=(1.0+eta2)/2.0;
	e2m=(1.0-eta2)/2.0;
	// somar sobre xi1 (i)
	aux=0.0;
	for(i=0;i<Q[0];i++){
	  // coordenadas
	  eta1=xGQ[0][i];
	  e1p=(1.0+eta1)/2.0;
	  e1m=(1.0-eta1)/2.0;
	  
	  x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
	  x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
	  x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
      
	  Fa=Psia(P[0], p,eta1);
	  
	  aux+=(wGQ[0][i]*Fa*func(x1,x2,x3)*JV[i+Q[0]*(j+Q[1]*k)]);// <==== Jacobian JV
	}
	fp[p][j][k]=aux;
      }
    }
  }
  // segunda matriz fpq(xi3)
  for(p=0;p<Pdim;p++){
    for(q=0;q<Pdim;q++){
      for(k=0;k<Q[2];k++){
	// somar sobre xi2 (j)
	aux=0.0;
	for(j=0;j<Q[1];j++){
	  eta2=xGQ[1][j];
	  Fb=Psib(P[0], P[1], p, q, eta2);
	  aux+=wGQ[1][j]*Fb*fp[p][j][k];
	}
	fpq[p][q][k]=aux;
      }
    }
  }
  for(n=0;n<nn;n++){
    p=mode_[n].p_val();
    q=mode_[n].q_val();
    r=mode_[n].r_val();
    aux=0.0;
    // fazer soma sobre xi3
    for(k=0;k<Q[2];k++){
      eta3=xGQ[2][k];
      Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
      aux+=wGQ[2][k]*Fc*fpq[p][q][k];
    }
    vec[n]=aux;
  }

  //for(n=0;n<nn;n++)vec[n]*=JV[0];//<========= Jacobiano e constante=====<
}
// *******************************************************************
// Mass Matrix entry M(m1, m2)
// *******************************************************************
// Evaluates the inner product of Phi_m1 by Phi_m2

double Tetrahedral::mass(int m1,int m2,const double JV[])
{
  int p1,q1,r1,p2,q2,r2;
  int i,j,k;
  double aux, aux1, aux2,Fa,Fb,Fc,Ga,Gb,Gc;
  //printf("Entrou mass(%d, %d)\n",m1,m2);

  // indexes of point m1
  p1=mode_[m1].p_val();
  q1=mode_[m1].q_val();
  r1=mode_[m1].r_val();
  // indexes of point m2
  p2=mode_[m2].p_val();
  q2=mode_[m2].q_val();
  r2=mode_[m2].r_val();
  aux=0.0;
  for(k=0;k<Q[2];k++){
    aux1=0.0;
    for(j=0;j<Q[1];j++){
      aux2=0.0;
      for(i=0;i<Q[0];i++){
	Fa=Psia(P[0], p1,xGQ[0][i]);
	Ga=Psia(P[0], p2,xGQ[0][i]);
	aux2+=wGQ[0][i]*Fa*Ga*JV[i+Q[0]*(j+Q[1]*k)];// <======Jacobiano===<
      }
      Fb=Psib(P[0], P[1], p1, q1, xGQ[1][j]);
      Gb=Psib(P[0], P[1], p2, q2, xGQ[1][j]);

      aux1+=wGQ[1][j]*Fb*Gb*aux2;
    }
    Fc=Psic(P[0], P[1], P[2], p1, q1, r1, xGQ[2][k]);
    Gc=Psic(P[0], P[1], P[2], p2, q2, r2, xGQ[2][k]);
 
    aux+=wGQ[2][k]*Fc*Gc*aux1;
  }
  return (aux);
};
// ************************************************************************
// Creates the local matrices using static condensation
// It uses NEWMAT matrices and operators
// Identical to Triangle::make_local_matrices()
// ************************************************************************
void Tetrahedral::make_local_matrices()
{
  int ni; 
  int i,ii,j,jj;
  ni=nn-nb;
  //printf("Making local matrices 1: Nb = %d  Ni= %d\n", Nb,Ni);
  int m1,m2;

  int p1,q1,r1,p2,q2,r2;
  int k;
  double aux, aux1, aux2,Fa,Fb,Fc,Ga,Gb,Gc;
  double M[nn][nn];
  for(m1=0;m1<nn;m1++){
    // indexes of point m1
    p1=mode_[m1].p_val();
    q1=mode_[m1].q_val();
    r1=mode_[m1].r_val();
    for(m2=m1;m2<nn;m2++){
      // indexes of point m2
      p2=mode_[m2].p_val();
      q2=mode_[m2].q_val();
      r2=mode_[m2].r_val();
      aux=0.0;
      for(k=0;k<Q[2];k++){
	aux1=0.0;
	for(j=0;j<Q[1];j++){
	  aux2=0.0;
	  for(i=0;i<Q[0];i++){
	    Fa=Psia(P[0], p1,xGQ[0][i]);
	    Ga=Psia(P[0], p2,xGQ[0][i]);
	    aux2+=wGQ[0][i]*Fa*Ga;
	  }
	  Fb=Psib(P[0], P[1], p1, q1, xGQ[1][j]);
	  Gb=Psib(P[0], P[1], p2, q2, xGQ[1][j]);
	  aux1+=wGQ[1][j]*Fb*Gb*aux2;
	}
	Fc=Psic(P[0], P[1], P[2], p1, q1, r1, xGQ[2][k]);
	Gc=Psic(P[0], P[1], P[2], p2, q2, r2, xGQ[2][k]);
	aux+=wGQ[2][k]*Fc*Gc*aux1;
      }
      M[m1][m2]=aux;
      M[m2][m1]=aux;
    }
  }

#ifdef _NEWMAT

  NEWMAT::Matrix mb(nb,nb);
  NEWMAT::Matrix mb_c(nb,nb);
  NEWMAT::Matrix mi(ni,ni);
  NEWMAT::Matrix mc(ni,nb);
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
//  printf("Saindo Tetrahedral::make_local_matrices: nb = %d  Ni= %d\n", nb,Ni);
#endif
};
 
/*
void Tetrahedral::make_mass(double ** MM, const double JV[])
{
  int ni; 
  int i,ii,j,jj;
  ni=nn-nb;
  int m1,m2;

  int p1,q1,r1,p2,q2,r2;
  int k;
  double aux, aux1, aux2,Fa,Fb,Fc,Ga,Gb,Gc;
  for(m1=0;m1<nn;m1++){
    // indexes of point m1
    p1=mode_[m1].p_val();
    q1=mode_[m1].q_val();
    r1=mode_[m1].r_val();
    for(m2=m1;m2<nn;m2++){
      // indexes of point m2
      p2=mode_[m2].p_val();
      q2=mode_[m2].q_val();
      r2=mode_[m2].r_val();
      aux=0.0;
      for(k=0;k<Q[2];k++){
	aux1=0.0;
	for(j=0;j<Q[1];j++){
	  aux2=0.0;
	  for(i=0;i<Q[0];i++){
	    Fa=Psia(P[0], p1,xGQ[0][i]);
	    Ga=Psia(P[0], p2,xGQ[0][i]);
	    aux2+=wGQ[0][i]*Fa*Ga*JV[i+Q[0]*(j+Q[1]*k)];//<===Jacobiano===<
	  }
	  Fb=Psib(P[0], P[1], p1, q1, xGQ[1][j]);
	  Gb=Psib(P[0], P[1], p2, q2, xGQ[1][j]);
	  aux1+=wGQ[1][j]*Fb*Gb*aux2;
	}
	Fc=Psic(P[0], P[1], P[2], p1, q1, r1, xGQ[2][k]);
	Gc=Psic(P[0], P[1], P[2], p2, q2, r2, xGQ[2][k]);
	aux+=wGQ[2][k]*Fc*Gc*aux1;
      }
      MM[m1][m2]=aux;
      MM[m2][m1]=aux;
    }
  }
}
*/
// ************************************************************************
/*
void Tetrahedral::make_mass_matrices(int NFields)
{
  // ***************************
  make_local_matrices();
  // ***************************
  if (NFields > 1)duplicar_mass_matrix(NFields);
#ifdef PRINTF_ON
  printf("Saindo Tetrahedral::make_mass_matrices(%2d): nb = %d  ni= %d\n",NFields,nb,nn-nb);
#endif
};
 */

/*
// ************************************************************************
// Reordena os nos e calcula f_mask do tetrahedro
// ************************************************************************
void Tetrahedral::ordenar4(int n[], int f_mask[])
{
// ********************************* Testada em 30/10/2016 *************
  int M1=0;
  int M2=0;
  int temp=n[3];
  int i;
  int n_old[4];
  int f_mask[4];
  for(i=0;i<4;i++)n_old[i]=n[i];
  for(i=0;i<4;i++){
    if(n[i]>M1) {
      M2=M1;
      M1=n[i];
    }
    else if(n[i]>M2) M2=n[i];
  }
  i=0;
  while(n[i]!=M1)i++;
  if(i!=3){
    if(i==0){
      n[3]=n[0];
      n[0]=n[1];
      n[1]=temp;
    }
    else if(i==1){
      n[3]=n[1];
      n[1]=n[2];
      n[2]=temp;
    }
    else {
      n[3]=n[2];
      n[2]=n[1];
      n[1]=temp;
    }
  }
  i=0;
  temp=n[2];
  while(n[i]!=M2)i++;
  if(i!=2){
    if(i==0){
      n[2]=n[0];
      n[0]=n[1];
      n[1]=temp;
    }
    else {
      n[2]=n[1];
      n[1]=n[0];
      n[0]=temp;
    }
  }
  // Renumerar as faces
  // Numeracao segue esquema do Gambit
  // face_numero = f_mask[face_original(fornecida pelo Gambit)]
  
  // Face 0 ABC recebe a face antiga
  if(n[3]==n_old[3])     f_mask[0]=0;
  else
    if(n[3]==n_old[2])   f_mask[1]=0;
    else
      if(n[3]==n_old[1]) f_mask[2]=0;
      else               f_mask[3]=0;
  
  // Face 1 ABD recebe a face antiga
  if(n[2]==n_old[3])     f_mask[0]=1;
  else
    if(n[2]==n_old[2])   f_mask[1]=1;
    else
      if(n[2]==n_old[1]) f_mask[2]=1;
      else               f_mask[3]=1;
  
  // Face 2 BCD recebe a face antiga
  if(n[1]==n_old[3])     f_mask[0]=2;
  else
    if(n[1]==n_old[2])   f_mask[1]=2;
    else
      if(n[1]==n_old[1]) f_mask[2]=2;
      else               f_mask[3]=2;
  
  // Face 3 ACD recebe a face antiga
  if(n[0]==n_old[3])     f_mask[0]=3;
  else
    if(n[0]==n_old[2])   f_mask[1]=3;
    else
      if(n[0]==n_old[1]) f_mask[2]=3;
      else               f_mask[3]=3;
// **********************************
 
};
*/

// ************************************************************************
// ************************************************************************
void Tetrahedral::printtofile(FILE * fout,const double u[], 
			      double (*func)(double,double,double),
			      const Vertice vert[], const int map[])
{
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  double e1p,e1m,e2p,e2m,e3p,e3m;
  double x1,x2,x3;//,y1,y2,y3,z1,z2,z3;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int m1;
  int i,j,k,p,q,r;
  //printf("Entrou em Tetrahedral printtofile\n");
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
  aux=0.0;
  
  for(k=0;k<Q[2];k++){
    eta3=xGQ[2][k];
    e3p=(1.0+eta3)/2.0;
    e3m=(1.0-eta3)/2.0;
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      e2p=(1.0+eta2)/2.0;
      e2m=(1.0-eta2)/2.0;
      for(i=0;i<Q[0];i++){
	// coordenadas
	eta1=xGQ[0][i];
	e1p=(1.0+eta1)/2.0;
	e1m=(1.0-eta1)/2.0;
	
	x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
	x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
	x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
	aux=0.0;
	for(m1=0;m1<nn;m1++){
	  p=mode_[m1].p_val();
	  q=mode_[m1].q_val();
	  r=mode_[m1].r_val();
	  Fa=Psia(P[0], p,eta1);
	  Fb=Psib(P[0], P[1], p, q, eta2);
	  Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
	  aux+=Fa*Fb*Fc*u[m1];;
	}
	fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
      }
    }
  }
};
// ************************************************************************
void Tetrahedral::printtofile(FILE * fout,const double u[], 
			      const Vertice vert[], const int map[])
{
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  double e1p,e1m,e2p,e2m,e3p,e3m;
  double x1,x2,x3;//,y1,y2,y3,z1,z2,z3;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int m1;
  int i,j,k,p,q,r;
  //printf("Entrou em Tetrahedral printtofile\n");
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
  aux=0.0;
  
  for(k=0;k<Q[2];k++){
    eta3=xGQ[2][k];
    e3p=(1.0+eta3)/2.0;
    e3m=(1.0-eta3)/2.0;
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      e2p=(1.0+eta2)/2.0;
      e2m=(1.0-eta2)/2.0;
      for(i=0;i<Q[0];i++){
	// coordenadas
	eta1=xGQ[0][i];
	e1p=(1.0+eta1)/2.0;
	e1m=(1.0-eta1)/2.0;
	
	x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
	x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
	x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
	aux=0.0;
	for(m1=0;m1<nn;m1++){
	  p=mode_[m1].p_val();
	  q=mode_[m1].q_val();
	  r=mode_[m1].r_val();
	  Fa=Psia(P[0], p,eta1);
	  Fb=Psib(P[0], P[1], p, q, eta2);
	  Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
	  aux+=Fa*Fb*Fc*u[m1];;
	}
	fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
      }
    }
  }
};
void Tetrahedral::printtoarray(const double u[], 
			       const Vertice vert[], const int map[],
			       double x[], double y[], double z[], 
			       double ftemp[])
{
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  double e1p,e1m,e2p,e2m,e3p,e3m;
 // double x1,x2,x3;//,y1,y2,y3,z1,z2,z3;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int m1;
  int i,j,k,p,q,r;
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
  aux=0.0;
  int count=0;
  for(k=0;k<Q[2];k++){
    eta3=xGQ[2][k];
    e3p=(1.0+eta3)/2.0;
    e3m=(1.0-eta3)/2.0;
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      e2p=(1.0+eta2)/2.0;
      e2m=(1.0-eta2)/2.0;
      for(i=0;i<Q[0];i++){
	// coordenadas
	eta1=xGQ[0][i];
	e1p=(1.0+eta1)/2.0;
	e1m=(1.0-eta1)/2.0;
	x[count]=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
	y[count]=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
	z[count]=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
	aux=0.0;
	for(m1=0;m1<nn;m1++){
	  p=mode_[m1].p_val();
	  q=mode_[m1].q_val();
	  r=mode_[m1].r_val();
	  Fa=Psia(P[0], p,eta1);
	  Fb=Psib(P[0], P[1], p, q, eta2);
	  Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
	  aux+=Fa*Fb*Fc*u[m1];;
	}
	ftemp[count++]=aux;
	//fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
      }
    }
  }
};

void Tetrahedral::evalGQ(double f0[],double f1[],
			 const double u0[],const double u1[])
{
  double aux0,aux1,Fa,Fb,Fc,eta1,eta2,eta3;
  int m1,n;
  int i,j,k,p,q,r;
  n=0;
  for(k=0;k<Q[2];k++){
    eta3=xGQ[2][k];
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      for(i=0;i<Q[0];i++){
	eta1=xGQ[0][i];
	aux0=0.0;
	aux1=0.0;
	for(m1=0;m1<nn;m1++){
	  p=mode_[m1].p_val();
	  q=mode_[m1].q_val();
	  r=mode_[m1].r_val();
	  Fa=Psia(P[0], p,eta1);
	  Fb=Psib(P[0], P[1], p, q, eta2);
	  Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
	  aux0+=Fa*Fb*Fc*u0[m1];
	  aux1+=Fa*Fb*Fc*u1[m1];
	}
	//fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
	f0[n]=aux0;
	f1[n++]=aux1;
      }
    }
  }
};

// Falta implementar
void Tetrahedral::evalGQ(double f0[],const double u0[],const int NF,const int nvar)
{

}
// ************************************************************************
// Evaluates the value of the field at the vertices of the element
// ************************************************************************
void Tetrahedral::computeVertice(double f_vert[],const double u[], 
				 const Vertice vert[], const int map[])
{
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  double e1p,e1m,e2p,e2m,e3p,e3m;
  double x1,x2,x3;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int m1;
  int i,p,q,r;
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

  i=0;
  for(eta3=-1.0;eta3<=1.0;eta3+=2.0){
    e3p=(1.0+eta3)/2.0;
    e3m=(1.0-eta3)/2.0;

    for(eta2=eta3;eta2<=1.0;eta2+=2.0){
      e2p=(1.0+eta2)/2.0;
      e2m=(1.0-eta2)/2.0;

      for(eta1=eta2;eta1<=1.0;eta1+=2.0){
	// coordenadas
	e1p=(1.0+eta1)/2.0;
	e1m=(1.0-eta1)/2.0;
	
	x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
	x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
	x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
	aux=0.0;
	for(m1=0;m1<nn;m1++){
	  p=mode_[m1].p_val();
	  q=mode_[m1].q_val();
	  r=mode_[m1].r_val();
	  Fa=Psia(P[0], p,eta1);
	  Fb=Psib(P[0], P[1], p, q, eta2);
	  Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
	  aux+=Fa*Fb*Fc*u[m1];;
	}
	//fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
	f_vert[map[i++]]=aux;
      }
    }
  }
};
// ************************************************************************
// Evaluates the value of the field at points
// ************************************************************************
void Tetrahedral::computeAtPoints(const int npoints, const double LocCoord[],const double u[], 
				  const Vertice vert[], const int map[],double f[],double GloCoord[])
{
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  double e1p,e1m,e2p,e2m,e3p,e3m;
  double x1,x2,x3;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int m1;
  int i,p,q,r;
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
    eta1=LocCoord[3*i  ];
    eta2=LocCoord[3*i+1];
    eta3=LocCoord[3*i+2];
    e3p=(1.0+eta3)/2.0;
    e3m=(1.0-eta3)/2.0;

    e2p=(1.0+eta2)/2.0;
    e2m=(1.0-eta2)/2.0;
    // coordenadas
    e1p=(1.0+eta1)/2.0;
    e1m=(1.0-eta1)/2.0;
    
    x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
    x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
    x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
    GloCoord[3*i  ]=x1;
    GloCoord[3*i+1]=x2;
    GloCoord[3*i+2]=x3;
    aux=0.0;
    for(m1=0;m1<nn;m1++){
      p=mode_[m1].p_val();
      q=mode_[m1].q_val();
      r=mode_[m1].r_val();
      Fa=Psia(P[0], p,eta1);
      Fb=Psib(P[0], P[1], p, q, eta2);
      Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
      aux+=Fa*Fb*Fc*u[m1];;
    }
    //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux, f(x1,x2,x3));
    f[i]=aux;
  }
};
// ************************************************************************
void Tetrahedral::Jacobian(const Vertice vert[],const int map[],double *JV)
// ************************************************************************
{
  double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  double xa,ya,xb,yb,xc,yc,xd,yd,za,zb,zc,zd;
  double Jacobian;
  // Calculo do Jacobiano do Tetraedro
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
  x1=xd-xa;
  x2=xb-xa;
  x3=xc-xa;
  y1=yd-ya;
  y2=yb-ya;
  y3=yc-ya;
  z1=zd-za;
  z2=zb-za;
  z3=zc-za;
  Jacobian=(x1*y2*z3 + y1*z2*x3 + z1*x2*y3 - x3*y2*z1 -y3*z2*x1 - z3*x2*y1)/8.0;
  if(Jacobian==0.0){
    printf("Jacobian=0\nVertices\n(%lf, %lf, %lf), (%lf, %lf, %lf), (%lf, %lf, %lf), (%lf, %lf, %lf)\n",xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd);
  }
#ifdef PRINTF_ON
  printf("J=%lf\n",Jacobian);
#endif
  for(int i=0;i<Q[0];i++)
    for(int j=0;j<Q[1];j++)
      for(int k=0;k<Q[2];k++)
	JV[i+Q[0]*(j+Q[1]*k)]=Jacobian;
};
// ************************************************************************

// ************************************************************************
// ************************************************************************

 // ************************************************************************
 
void Tetrahedral::teste_face(int a,int b,int c, int& sign, int& ng,
			     int& NG,int& NF,int Face[],int Fng[])
{
  int p=P[0];
  int ii,NF3=0;
  int flag; 
  
  sign=1;
  
  // Verificar se a face ja foi numerada
  flag=0; // reset flag; flag = 0 equivale a face nao numerada
  for (int i=0; i<NF && flag==0; i++){
    ii=3*i;
    if(Face[ii+2]==c){// Face ja numerada
      if(Face[ii]==a && Face[ii+1]==b){// Face percorrida no mesmo sentido
	ng= Fng[i];
	flag=1;
      }
      else if(Face[ii]==b && Face[ii+1]==a){//  face ja percorrida em sentido oposto
	ng= Fng[i];
	flag=1;
	sign=-1;
      }
    }
  }
  if(flag==0){// Face ainda nao numerada
    ng=NG;
    NF3=3*NF;
    Face[NF3]=a;
    Face[NF3+1]=b;
    Face[NF3+2]=c;
    Fng[NF]=ng;
    NF++; // acrescenta o numero de faces armezandas
    NG+=(p-1)*(p-2)/2;
  }
  if(NF3 > MAXNL || NG > MAXNG){
    printf("Ajuste indices\nNF3 = %d (MAXNL = %d)\nNG = %d (MAXNG = %d)\n",NF3,MAXNL,NG,MAXNG);
  }
};
// ************************************************************************
void Tetrahedral::map_aresta(int imin,int imax,int sign0,int ng0,int map[],int sgn[])
{
  int i,temp;
  if(sign0==-1){
    temp=1;
    for(i=imin; i<=imax;i++){
      map[i]=ng0;
      sgn[i]=temp;
      temp*=-1;// muda o sinal
      ng0++;
    }
  }
  else 
    for(i=imin; i<=imax;i++){
      map[i]=ng0;
      sgn[i]=1;
      ng0++;
    }
};
// ************************************************************************
void Tetrahedral::map_face(int & imin,int P0, int P1,
			   int sign0,int ng0,int map[],int sgn[])
{
  int a=imin; 
  int temp=1;
  for(int i=1; i<P0-1; i++){
    for(int j=1; j<P1-i;j++){
      map[a]=ng0;
      sgn[a]=temp;
      ng0++;
      a++;
    }
    temp*=sign0;// muda o sinal se sign0==-1
  }
  imin=a;// returns the next available address
};


// ************************************************************************
// Calculates the gradient wrt x1, x2, x3
// grad[m][i] = component i of the gradient at the point m(pqr)
// ************************************************************************
// Input function in vector form fvec
// ************************************************************************
void Tetrahedral::Gradiente(FILE * fout, double ** grad,  
			    const double fvec[],
			    const Vertice vert[], 
			    const int map[])
{
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,eta1,eta2,eta3,x1,x2,x3;
  double e3p, e3m,e2p,e2m,e1p,e1m;
  int i,j,k,l,m;
  double df[MAXQ*MAXQ*MAXQ][3];
  // double gvec[MAXQ*MAXQ*MAXQ];
  double a11,a12,a13,a21,a22,a23,a31,a32,a33, J3D;
  double b[3][3];
  double aux,aux0,aux1,aux2;
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
  // calculo dos coeficientes aij
  a11=0.5*(xb-xa);
  a12=0.5*(xc-xa);
  a13=0.5*(xd-xa);
  a21=0.5*(yb-ya);
  a22=0.5*(yc-ya);
  a23=0.5*(yd-ya);
  a31=0.5*(zb-za);
  a32=0.5*(zc-za);
  a33=0.5*(zd-za);
  // calculo de J3D
  J3D=a11*(a22*a33-a23*a32) - a12*(a21*a33-a23*a31) + a13*(a21*a32-a22*a31);
  // calculo da matriz b[i][j]
  b[0][0]=(a22*a33-a23*a32)/J3D;
  b[0][1]=-(a12*a33-a13*a32)/J3D;
  b[0][2]=(a12*a23-a13*a22)/J3D;
  b[1][0]=-(a21*a33-a23*a31)/J3D;
  b[1][1]=(a11*a33-a13*a31)/J3D;
  b[1][2]=-(a11*a23-a13*a21)/J3D;
  b[2][0]=(a21*a32-a22*a31)/J3D;
  b[2][1]=-(a11*a32-a12*a31)/J3D;
  b[2][2]=(a11*a22-a12*a21)/J3D;
  // calculo das derivadas com relacao a eta1, eta2, eta3
  for(int r=0;r<Q[2];r++){
    for(int q=0;q<Q[1];q++){
      for(int p=0;p<Q[0];p++){
	m=p+q*Q[0]+r*Q[0]*Q[1];
	aux0=0.0;
	aux1=0.0;
	aux2=0.0;
	for(l=0;l<Q[0];l++){
	  aux0+=D[p][l][0]*fvec[(l+q*Q[0]+r*Q[0]*Q[1])];
	  aux1+=D[q][l][1]*fvec[(p+l*Q[0]+r*Q[0]*Q[1])];
	  aux2+=D[r][l][2]*fvec[(p+q*Q[0]+l*Q[0]*Q[1])];
	}
	df[m][0]=aux0;
	df[m][1]=aux1;
	df[m][2]=aux2;
      }
    }
  }
  // calculo das derivadas com relacao a xi1, xi2, xi3
  for(int r=0;r<Q[2];r++){
    eta3=xGQ[2][r];
    for(int q=0;q<Q[1];q++){
      eta2=xGQ[1][q];
      for(int p=0;p<Q[0];p++){
	eta1=xGQ[0][p];
	m=p+q*Q[0]+r*Q[0]*Q[1];
	aux=(1-eta2)*(1-eta3);
  	aux0=(4/aux)*df[m][0]; // d/dxi_1
	aux1=2*(1+eta1)/aux*df[m][0] + 2/(1-eta3)*df[m][1];// d/dxi_2
	aux2=2*(1+eta1)/aux*df[m][0] + (1+eta2)/(1-eta3)*df[m][1] + df[m][2];// d/dxi_3
	df[m][0]=aux0;// d/dxi_1
	df[m][1]=aux1;// d/dxi_2
	df[m][2]=aux2;// d/dxi_3
      }
    }
  }
  // calculo das derivadas com relacao a x1, x2, x3
  for(int r=0;r<Q[2];r++){
    for(int q=0;q<Q[1];q++){
      for(int p=0;p<Q[0];p++){
	m=p+q*Q[0]+r*Q[0]*Q[1];
	for(i=0;i<3;i++){
	  aux=0;
	  for(j=0;j<3;j++){
	    aux+=b[j][i]*df[m][j];
	  }
	  grad[i][m]=aux;
	}
      }
    }
  }
   // Cheque do gradiente
   for(k=0;k<Q[2];k++){
     eta3=xGQ[2][k];
     e3p=(1.0+eta3)/2.0;
     e3m=(1.0-eta3)/2.0;
     for(j=0;j<Q[1];j++){
       eta2=xGQ[1][j];
       e2p=(1.0+eta2)/2.0;
       e2m=(1.0-eta2)/2.0;
       for(i=0;i<Q[0];i++){
	 eta1=xGQ[0][i];
 	e1p=(1.0+eta1)/2.0;
 	e1m=(1.0-eta1)/2.0;
 	// coordenadas x1, x2, x3
 	x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
 	x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
 	x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
 	
 	m=i+j*Q[0]+k*Q[0]*Q[1];
 	fprintf(fout,"testex %11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,g(x1,x2,x3),grad[0][m]);// Derivada com relacao a x no ponto de integracao m=m(i,j,k).
       }
     }
   }

};
// ************************************************************************
void Tetrahedral::Gradiente(double ** grad,  
			    const double fvec[], 
			    const Vertice vert[], 
			    const int map[])
{
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,eta1,eta2,eta3;
 
  int i,j,l,m;
  double df[MAXQ*MAXQ*MAXQ][3];
  
  double a11,a12,a13,a21,a22,a23,a31,a32,a33, J3D;
  double b[3][3];
  double aux,aux0,aux1,aux2;
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
  // calculo dos coeficientes aij
  a11=0.5*(xb-xa);
  a12=0.5*(xc-xa);
  a13=0.5*(xd-xa);
  a21=0.5*(yb-ya);
  a22=0.5*(yc-ya);
  a23=0.5*(yd-ya);
  a31=0.5*(zb-za);
  a32=0.5*(zc-za);
  a33=0.5*(zd-za);
  // calculo de J3D
  J3D=a11*(a22*a33-a23*a32) - a12*(a21*a33-a23*a31) + a13*(a21*a32-a22*a31);
  // calculo da matriz b[i][j]
  b[0][0]=(a22*a33-a23*a32)/J3D;
  b[0][1]=-(a12*a33-a13*a32)/J3D;
  b[0][2]=(a12*a23-a13*a22)/J3D;
  b[1][0]=-(a21*a33-a23*a31)/J3D;
  b[1][1]=(a11*a33-a13*a31)/J3D;
  b[1][2]=-(a11*a23-a13*a21)/J3D;
  b[2][0]=(a21*a32-a22*a31)/J3D;
  b[2][1]=-(a11*a32-a12*a31)/J3D;
  b[2][2]=(a11*a22-a12*a21)/J3D;
  // calculo das derivadas com relacao a eta1, eta2, eta3
  for(int r=0;r<Q[2];r++){
    for(int q=0;q<Q[1];q++){
      for(int p=0;p<Q[0];p++){
	m=p+q*Q[0]+r*Q[0]*Q[1];
	aux0=0.0;
	aux1=0.0;
	aux2=0.0;
	for(l=0;l<Q[0];l++){
	  aux0+=D[p][l][0]*fvec[(l+q*Q[0]+r*Q[0]*Q[1])];
	  aux1+=D[q][l][1]*fvec[(p+l*Q[0]+r*Q[0]*Q[1])];
	  aux2+=D[r][l][2]*fvec[(p+q*Q[0]+l*Q[0]*Q[1])];
	}
	df[m][0]=aux0;
	df[m][1]=aux1;
	df[m][2]=aux2;
      }
    }
  }
  // calculo das derivadas com relacao a xi1, xi2, xi3
  for(int r=0;r<Q[2];r++){
    eta3=xGQ[2][r];
    for(int q=0;q<Q[1];q++){
      eta2=xGQ[1][q];
      for(int p=0;p<Q[0];p++){
	eta1=xGQ[0][p];
	m=p+q*Q[0]+r*Q[0]*Q[1];
	aux=(1-eta2)*(1-eta3);
  	aux0=(4/aux)*df[m][0]; // d/dxi_1
	aux1=2*(1+eta1)/aux*df[m][0] + 2/(1-eta3)*df[m][1];// d/dxi_2
	aux2=2*(1+eta1)/aux*df[m][0] + (1+eta2)/(1-eta3)*df[m][1] + df[m][2];// d/dxi_3
	df[m][0]=aux0;// d/dxi_1
	df[m][1]=aux1;// d/dxi_2
	df[m][2]=aux2;// d/dxi_3
      }
    }
  }
  // calculo das derivadas com relacao a x1, x2, x3
  for(int r=0;r<Q[2];r++){
    for(int q=0;q<Q[1];q++){
      for(int p=0;p<Q[0];p++){
	m=p+q*Q[0]+r*Q[0]*Q[1];
	for(i=0;i<3;i++){
	  aux=0;
	  for(j=0;j<3;j++){
	    aux+=b[j][i]*df[m][j];
	  }
	  grad[i][m]=aux;
	}
      }
    }
  }
};
// ************************************************************************
// Calculates the gradient wrt x1, x2, x3
// grad[m][i] = component i of the gradient at the point m(pqr)
// ************************************************************************
// Input function in func(double,double,double)
// ************************************************************************
void Tetrahedral::Gradiente(FILE * fout,
			    double ** grad,
			    double (*func)(double, double, double),
			    const Vertice vert[],
			    const int map[])
{
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,eta1,eta2,eta3,x1,x2,x3;
  double e3p, e3m,e2p,e2m,e1p,e1m;
  int i,j,k,m;
  double fvec[MAXQ*MAXQ*MAXQ];
  
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
  // calculo do vetor contendo a funcao nos pontos de integracao de Gauss
  for(k=0;k<Q[2];k++){
    eta3=xGQ[2][k];
    e3p=(1.0+eta3)/2.0;
    e3m=(1.0-eta3)/2.0;
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      e2p=(1.0+eta2)/2.0;
      e2m=(1.0-eta2)/2.0;
      for(i=0;i<Q[0];i++){
	eta1=xGQ[0][i];
	e1p=(1.0+eta1)/2.0;
	e1m=(1.0-eta1)/2.0;
	// coordenadas x1, x2, x3
	x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
	x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
	x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
	
	m=i+j*Q[0]+k*Q[0]*Q[1];
	fvec[m]=func(x1,x2,x3);
      }
    }
  }
  Gradiente(fout, grad, fvec, vert, map); 
};
      
void Tetrahedral::print_nome(FILE * fout)
{
  fprintf(fout,"ELEMENTO TETRAEDRICO\n");
};
// ************************************************************************
// Mapear a numeracao do triangulo na face fnum do tetraedro sobre a numeracao
// do triangulo padrao que e retornada em trimap[] e marcar o no global em bflag com
// valor 0 (conhecido)
// ************************************************************************
void Tetrahedral::localFaceModeMap(const int fnum,
                                   const int P,
                                   int trimap[])
{
  int mar[6][P-1];
  int mface[4][(P-1)*(P-2)/2];
  int ar[3]; // aresta que compoem a face
  
  // criar numeracao local dos modos do tetraedro
  int a=4; // numero do proximo modo
  for(int i=0;i<6;++i) { // loop sobre as arestas
    for(int j=0;j<P-1;++j)
      mar[i][j] = a++;
  }
  for(int k=0;k<4;++k) { // loop sobre as faces
    int count=0;
    for(int i=1;i<P;++i)
      for(int j=1;j<P-i;++j)
        mface[k][count++] = a++;
  }
  
  int l=0;
 
  // Flag dos nos
  for(int i=0;i<3;++i) {
    trimap[l++] = face[fnum][i];
    ar[i] = face_aresta[fnum][i];
  }
  for(int i=0;i<3; ++i) { //loop sobre as arestas do triangulo
    for(int j=0;j<P-1;++j)
      trimap[l++] = mar[ar[i]][j];
  }
  int aux = (P-1)*(P-2)/2;
  for(int j=0;j<aux;++j)
    trimap[l++] = mface[fnum][j];

#ifdef PRINTF_ON
  printf("trimap da face %d\n",fnum);
  for(int i=0;i<l;++i){
    printf("trimap[%d] = %d ", i,trimap[i]);
    mode_[trimap[i]].print();
  }

 printf("Saindo de Tetrahedral::set: nv = %d  nn = %d (< MAXMODES = %d) nb = %d (< MAXNB = %d)\n",nv,nn,MAXMODES,nb,MAXNB);
#endif
};
// ************************************************************************
// ************************************************************************
void Tetrahedral::Dirichlet(const int face_num,
                            const Vertice vert[],
                            const int vert_map[], // numero global dos vertices dos nos
                            const int nmap[], // mapeamento global dos modos locais
                            const int sgn[],
                            int bflag[],
                            double Xbc[],
                            double (*func)(double,double,double))
{
  
   //cout << "Entrando em Tetrahedral::Dirichlet para a face "<< face_num << endl;
 
  const int varn = 0;
  //const int NFields = 1;
 // double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
 // Vertice lvert[3]; //Vertices para problema local
  int lvert_map[3]={0,1,2};
  int trimap[MAXMODES];// mapeia modos do triangulo sobre os do tetraedro
  //double Xbc1[MAXNB];
  //int bflag1[MAXNB];
  int p0,p1;
  int i,ii;
   // construcao do mapa dos vertices locais em vertices globais
  for(int i=0;i<3;++i) {
    lvert_map[i]=vert_map[face[face_num][i]];
  }
  p0=P[fd0[face_num]];
  p1=P[fd1[face_num]];
  
  int q=Q[0];
  //dados do triangulo (stdel)
  Triangle * triang = new Triangle(p0,q);
  PhElem<1> * localphel = new PhElem<1>();// So um campo
  localphel->set_ptr_stdel(triang,triang);
  localphel->set_type(2); // triangulo
  localphel->set_ptvert(vert); // array de vertices recebido nos argumentos
  localphel->set_Vert_map(3,lvert_map);
  //localphel->compute_JV(0);// Jacobiano J deve ser calculado antes do vetor
  localphel->inicia_vetores();// zera o vetor b de localphel
  int count = 0;
  localphel->inicia_gbnmap(0,count);// inicia gbnmap[i]=i,sgn[i]=1, i=0,nn-1;
  // Encontrar os nos.
  //Criar trimap[]; mapeamento dos nos do triangulo standard sobre 
  // os nos globais 
  localFaceModeMap(face_num,p0,trimap); // mapeia o triangulo local no tetraedro
 
  //cout << "chamando projetar_C0 para triangulo (localphel)"<< endl;
  localphel->projetar_C0(nullptr, func, varn);
 
   //mapeia os resultados sobre os modos globais
  for(i=0;i<triang->nn_val();++i){
    ii=trimap[i]; // indice do modo no tetraedro
    Xbc[nmap[ii]]=localphel->show_u(varn,i)*sgn[ii];
    bflag[nmap[ii]]=0; // conhecido
  }
  localphel->finaliza_vetores();
  delete localphel;
  delete triang;
  //printf("Saindo Tetrahedral::Dirichlet para face %d\n\n",face_num);
};

// *****************************************************************
void Tetrahedral::eval_Phi(const int n,double Phi[])
{
  //printf("Hexahedral::eval_Phi\n");
  int q2=Q[0]*Q[1]*Q[2];
  for(int i=0;i<q2;++i)Phi[i]=Phi_val[n][i];
};
// void  Tetrahedral::Jacobian_Vector(const Vertice vert[],const int map[],double JV[])
// {
//   int q2=Q[0]*Q[1]*Q[2];
//   double J=Jacobian(vert,map);
//   for(int i=0;i<q2;i++)JV[i]=J;
// };

// Sem implementacao

void Tetrahedral::elem_traces(const Vertice vert[],const int map[],const int sinal[],
                           double *** TP,double **** TGP,double * Jb)
{};
void Tetrahedral::eval_GradPhi(const Vertice vert[], const int map[],const int m,double ** der)
{
  printf("Nao implementado\n");
};

void Tetrahedral::computeFuncGQ(double f_[],
                             const Vertice vert[], const int map[],
                             double (*func)(double,double,double))
{
  double eta1,eta2,eta3;
  double e1p,e1m,e2p,e2m,e3p,e3m;
  double x1,x2,x3;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int i,j,k;
  //printf("Entrou em Tetrahedral printtofile\n");
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
  int count=0;

  for(k=0;k<Q[2];k++){
    eta3=xGQ[2][k];
    e3p=(1.0+eta3)/2.0;
    e3m=(1.0-eta3)/2.0;
    for(j=0;j<Q[1];j++){
      eta2=xGQ[1][j];
      e2p=(1.0+eta2)/2.0;
      e2m=(1.0-eta2)/2.0;
      for(i=0;i<Q[0];i++){
      // coordenadas
        eta1=xGQ[0][i];
        e1p=(1.0+eta1)/2.0;
        e1m=(1.0-eta1)/2.0;
              
        x1=((xa*e1m+xb*e1p)*e2m+e2p*xc)*e3m+e3p*xd;
        x2=((ya*e1m+yb*e1p)*e2m+e2p*yc)*e3m+e3p*yd;
        x3=((za*e1m+zb*e1p)*e2m+e2p*zc)*e3m+e3p*zd;
        
        f_[count++]=func(x1,x2,x3);
      }
    }
  }
};

void Tetrahedral::printGQtofile(FILE * fout,const double ftemp[],
                             const double ftemp1[],
                             const Vertice vert[], const int map[])
{};

void Tetrahedral::printwGQtofile(FILE * fout,
                              const Vertice vert[],
                              const int map[],
                              const double JV[])
{};
void Tetrahedral::make_mass_matrices(int NFields)
{/*
  //printf("Tetrahedral make_mass_matrices\n");
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

void Tetrahedral::vector_of_integral_of_f_Phi_dv(double vec[],
                                    const double func[],
                                    const double JV[])
{
  int n,p,q,r;
  int i,j,k;
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  int Pdim=P[0]+2;
  double fp[Pdim][Q[1]][Q[2]], fpq[Pdim][Pdim][Q[2]];
  
  // primeira matriz fp(xi2,xi3)
  for(p=0;p<Pdim;p++){
    // fazer loop sobre xi2 e xi3
    for(k=0;k<Q[2];k++){
   
      for(j=0;j<Q[1];j++){
     
        // somar sobre xi1 (i)
        aux=0.0;
        for(i=0;i<Q[0];i++){
          // coordenadas
          eta1=xGQ[0][i];
          
          Fa=Psia(P[0], p,eta1);
          n = i + Q[0]*(j + Q[1]*k);
          aux+=(wGQ[0][i]*Fa*func[n]*JV[n]);// <==== Jacobian JV
        }
        fp[p][j][k]=aux;
      }
    }
  }
  // segunda matriz fpq(xi3)
  for(p=0;p<Pdim;p++){
    for(q=0;q<Pdim;q++){
      for(k=0;k<Q[2];k++){
        // somar sobre xi2 (j)
        aux=0.0;
        for(j=0;j<Q[1];j++){
          eta2=xGQ[1][j];
          Fb=Psib(P[0], P[1], p, q, eta2);
          aux+=wGQ[1][j]*Fb*fp[p][j][k];
        }
        fpq[p][q][k]=aux;
      }
    }
  }
  for(n=0;n<nn;n++){
    p=mode_[n].p_val();
    q=mode_[n].q_val();
    r=mode_[n].r_val();
    aux=0.0;
    // fazer soma sobre xi3
    for(k=0;k<Q[2];k++){
      eta3=xGQ[2][k];
      Fc=Psic(P[0],P[1],P[2],p,q,r,eta3);
      aux+=wGQ[2][k]*Fc*fpq[p][q][k];
    }
    vec[n]=aux;
  }
};

void Tetrahedral::teste(int & v)
{
  v=100000;
};

void Tetrahedral::trace(const int lado, const int qmax, const int sinal,
                     const double * valores, // valores nos pontos de Gauss
                     double * saida)
{};
const int Tetrahedral::aresta_lvert(const int & i, const int & j) const {return aresta[i][j];};
const int Tetrahedral::face_lvert(const int & i, const int & j) const {return face[i][j];};
const int Tetrahedral::show_nvf(const int &i) const {return nvf[i];};
const int Tetrahedral::show_face_tipo(const int &i) const {return face_tipo[i];};
const int Tetrahedral::show_fd0(const int &i) const {return fd0[i];};
const int Tetrahedral::show_fd1(const int &i) const {return fd1[i];};
const int Tetrahedral::show_fv2(const int &i) const {return fv2[i];};
const int Tetrahedral::show_ind_mode(const int & i, const int & j, const int & k ) const {return ind_mode_[i][j][k];};
void Tetrahedral::superficie_externa(const int Vert_map[],const Vertice vert[],
                                     const int & num_local,
                                     double & area,double normal[3])
{
  int v0,v1,v2;
  v0=Vert_map[face[num_local][0]];
  v1=Vert_map[face[num_local][1]];
  v2=Vert_map[face[num_local][2]];
  double l0x,l0y,l0z;
  double l1x,l1y,l1z;
  l0x=vert[v1].x - vert[v0].x;
  l0y=vert[v1].y - vert[v0].y;
  l0z=vert[v1].z - vert[v0].z;
  
  l1x=vert[v2].x - vert[v0].x;
  l1y=vert[v2].y - vert[v0].y;
  l1z=vert[v2].z - vert[v0].z;

  // Produto vetorial das arestas vezes o sinal
  normal[0] = sinal_normal[num_local] * (l0y*l1z - l1y*l0z);
  normal[1] = sinal_normal[num_local] * (l0z*l1x - l1z*l0x);
  normal[2] = sinal_normal[num_local] * (l0x*l1y - l1x*l0y);
  
  area = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  
  normal[0]/=area;
  normal[1]/=area;
  normal[2]/=area;
};
