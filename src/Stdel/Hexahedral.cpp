#include "spectral.h"
//// #include "PhElem.cpp"
//// #include "PhElem_locais.cpp"
// ************************************************************************
// Class Hexahedral
// ************************************************************************
// ************************************************************************
Hexahedral::Hexahedral(int p0,int q0)
{
  //printf("Hexahedral::Hexahedral\n");
  tipo=5;
  ndim=3;
  nv=8;
  ne=12;
  nf=6;
  nborder=6;
  vtk_type=12;
  emapi = new int [ne+1];
  emapv = new int [ne*(p0+1)];
  set(p0, q0);
};
// ************************************************************************
Hexahedral::~Hexahedral()
{
  //printf("Hexahedral::~Hexahedral\n");
  //printf("Destruindo Hexahedral\n");
  delete [] emapi; emapi=nullptr;
  delete [] emapv; emapv=nullptr;
  delete [] D_Phi_val; D_Phi_val=nullptr;
  
  //libera memoria de ind_mode_
  for (int i=0;i<=P[0];++i){
    for(int j=0;j<=P[1];++j){
      delete [] ind_mode_[i][j];  ind_mode_[i][j] = nullptr;
    }
    delete [] ind_mode_[i];  ind_mode_[i] = nullptr;
  }
  delete [] ind_mode_;  ind_mode_ = nullptr;
};
// ************************************************************************

// Inicializa os dados do elemento
void Hexahedral::set(int p0, int q0)
{
  //printf("Hexahedral::set\n");
 // int i,j,k;
  // Mapeamento das arestas
  // Indices que marcam o inicio das arestas
  for (int i=0; i <= 12; ++i) {
    emapi[i]=i*(p0+1);
  }

  P[0]=p0; Q[0]=q0; gqt[0]=3;// Gauss-Lobatto-Jacobi
  P[1]=p0; Q[1]=q0; gqt[1]=3;// Gauss-Lobatto-Jacobi
  P[2]=p0; Q[2]=q0; gqt[2]=3;// Gauss-Lobatto-Jacobi
  NGQP=Q[0]*Q[1]*Q[2];
  qborder = q0*q0;
    // incremento no indice dos pontos de Gauss na direcao i
    inc[0]=1;
    inc[1]=Q[0];
    inc[2]=Q[0]*Q[1];

  //aloca memoria para ind_mode_
  ind_mode_ = new int ** [P[0]+1];
  for (int i=0;i<=P[0];++i){
    ind_mode_[i] = new int * [P[1]+1];
    for(int j=0;j<=P[1];++j){
      ind_mode_[i][j] = new int [P[2]+1];
    }
  }

#ifdef PRINTF_ON
  printf("Em Hexahedral::set(%d,%d)\n",p0,q0);
  printf("Em Hexahedral::set apos Sair de Hexahedral::setStdel\n");
  printf("P: %d %d\nQ: %d %d\ngqt : %d %d\n",P[0],P[1],Q[0],Q[1],gqt[0],gqt[1]);
#endif
  // Ordem local

  //     C ---- D
  //   / |     /|
  // G ------H  |
  // |   A --|- B
  // | /     | /
  // E-------F
  // Modos dos vertices
  mode_[0].set_mode(0   ,0   ,0);      // A
  ind_mode_[0][0][0] = 0;

  mode_[1].set_mode(P[0],0   ,0);      // B
  ind_mode_[P[0]][0][0] = 1;

  mode_[2].set_mode(0   ,P[1],0);      // C
  ind_mode_[0][P[1]][0] = 2;

  mode_[3].set_mode(P[0],P[1],0);      // D
  ind_mode_[P[0]][P[1]][0] = 3;

  mode_[4].set_mode(0   ,0   ,P[2]);   // E
  ind_mode_[0][0][P[2]] = 4;

  mode_[5].set_mode(P[0],0   ,P[2]);   // F
 ind_mode_[P[0]][0][P[2]] = 5;

  mode_[6].set_mode(0   ,P[1],P[2]);   // G
  ind_mode_[0][P[1]][P[2]] = 6;

  mode_[7].set_mode(P[0],P[1],P[2]);   // H
  ind_mode_[P[0]][P[1]][P[2]] = 7;
  int a=8;

  int dir[3], v[2];
  int vv[3] = {1,2,4};
  int ind[3];
  
  for(int i=0;i<12;++i) { // loop sobre as arestas
    dir[2] = (aresta[i][1]-aresta[i][0])/2;
    dir[0] = (dir[2] + 1) % 3;
    dir[1] = (dir[2] + 2) % 3;
    v[0] = (vv[dir[0]] & aresta[i][0])/vv[dir[0]];
    v[1] = (vv[dir[1]] & aresta[i][0])/vv[dir[1]];
   
    ind[dir[0]] = (v[0] == 1) ? P[dir[0]] : 0;
    ind[dir[1]] = (v[1] == 1) ? P[dir[1]] : 0;

    int k=emapi[i]; // aresta i
    emapv[k]=aresta[i][0]; // ponto inicial
    k++;
    for(int p=1;p < P[dir[2]];++p) {
      ind[dir[2]]= p;
      mode_[a].set_mode(ind[0],ind[1],ind[2]);
      ind_mode_[ind[0]][ind[1]][ind[2]] = a;
      emapv[k]=a;
      k++;
      a++;
    }
    emapv[k]=aresta[i][1]; // ponto final
  }

  //fazer o loop sobre as faces
  for(int i=0; i < nf; ++i){
    // int sgn[2]={1,1}; // os sinais no elemento padrao sao sempre positivos
    
    // determina as direcoes
    dir[0]=(face[i][1]-face[i][0]);
    dir[1]=(face[i][3]-face[i][0]);
    
    // os sinais no elemento padrao sao sempre positivos
    //if(dir[0] < 0) {sgn[0] = -1;  dir[0] *= -1;}
    //if(dir[1] < 0) {sgn[1] = -1;  dir[1] *= -1;}
    
    dir[2]=((dir[0] ^ 7) ^ dir[1]);
    
    int v2=(face[i][0] & face[i][1] & face[i][2] & face[i][3] & dir[2])/dir[2];
    
    dir[0] /=2;
    dir[1] /=2;
    dir[2] /=2;

    ind[dir[2]] = (v2 == 1) ? P[dir[2]] : 0;
    for(int p=1;p < P[dir[0]];++p) {
      ind[dir[0]]= p;
      for(int q=1;q < P[dir[1]];++q) {
        ind[dir[1]]= q;
        mode_[a].set_mode(ind[0],ind[1],ind[2]);
        ind_mode_[ ind[0] ][ ind[1] ][ ind[2] ] = a;
        ++a;
      }
    }
  }
  //
  nb=a;
  
  //Verificar
  // Interior modes : j runs fastest
  for(int i=1; i<P[0]; ++i){
    for(int j=1; j<P[1]; ++j){
      for(int k=1; k<P[2]; ++k) {
          mode_[a].set_mode(i,j,k);
          ind_mode_[ i ][ j ][ k ] = a;
          ++a;
      }
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
  printf("Em Hexahedral::set  passou gauss_parameters\n");
 // printf("nb=%d emapv tem %d valores\n",nb,k);
  printf("Saindo de Hexahedral::set: nv = %d  nn = %d (< MAXMODES = %d) nb = %d (< MAXNB = %d)\n",nv,nn,MAXMODES,nb,MAXNB);
  printf("Hexahedral::set nn = %d nb= %d\n",nn,nb);
#endif
 
  // Construcao da matriz Phi_val[nn][Q[0]*Q[1]*Q[2]]
  double eta1,eta2,eta3;
  double Fa[Q[0]],Fb[Q[1]],Fc[Q[2]];
  for(int m=0;m<nn;++m){
    
    int p= mode_[m].p_val();
    int q= mode_[m].q_val();
    int r= mode_[m].r_val();
    for(int k=0;k<Q[2];++k){
      eta3=xGQ[2][k];
      Fc[k]=Psia(P[2],r,eta3);
    }
    for(int j=0;j<Q[1];++j){
      eta2=xGQ[1][j];
      Fb[j]=Psia(P[1],q,eta2);
    } 
    for(int i=0;i<Q[0];++i){
      eta1=xGQ[0][i];
      Fa[i]=Psia(P[0],p,eta1);
    }
    
    int n=0;
    for(int k=0;k<Q[2];++k){
      for(int j=0;j<Q[1];++j){
        for(int i=0;i<Q[0];++i){// i runs fastest
          Phi_val[m][n++]=Fa[i]*Fb[j]*Fc[k];
        }
      }
    }
  }
 FILE * fout = fopen("hexa.txt", "wb");
 print(fout);
 fclose(fout);
};
// 01/07/2014
// ************************************************************************
void Hexahedral::print(FILE * fout)
{
  printStdel();
	//printf("Hexahedral::print\n");
  fprintf(fout,"Hexahedral::print\n");
  for(int i=0;i<nn; i++){
    fprintf(fout,"no= %4d ",i);
    mode_[i].print(fout);
    for(int k =0; k<NGQP_val();++k)
      fprintf(fout,"nGQ = %4d Phi_val = %g\n",k,Phi_val[i][k]);
  }
  
  //fprintf(fout,"\n");
};

// ************************************************************************
double Hexahedral::mass(int m1,int m2,const double JV[])
{
  //printf("Hexahedral::mass\n");
  int p1,q1,r1,p2,q2,r2;
  //int flag1, flag2;
  int i,j,k;
  double aux, aux0,aux1, Fa,Fb,Ga,Gb,Fc,Gc;
  // printf("Entrou Hexahedral::mass(%d, %d, JV[])\n",m1,m2);
  // indexes of point m1
  p1= mode_[m1].p_val();
  q1= mode_[m1].q_val();
  r1= mode_[m1].r_val();
  // indexes of point m2
  p2= mode_[m2].p_val();
  q2= mode_[m2].q_val();
  r2= mode_[m2].r_val();
  aux=0.0;
  for(k=0;k<Q[2];++k){
    aux1=0.0;
    for(j=0;j<Q[1];++j){
      aux0=0.0;
      for(i=0;i<Q[0];++i){
        Fa=Psia(P[0], p1,xGQ[0][i]);//primeiro
        Ga=Psia(P[0], p2,xGQ[0][i]);
        aux0+=wGQ[0][i]*Fa*Ga*JV[i+Q[0]*(j+Q[1]*k)];
      }
      Fb=Psia(P[1], q1, xGQ[1][j]);
      Gb=Psia(P[1], q2, xGQ[1][j]);
      aux1+=wGQ[1][j]*Fb*Gb*aux0;
    }
    Fc=Psia(P[2], r1, xGQ[2][k]);
    Gc=Psia(P[2], r2, xGQ[2][k]);
    aux+=wGQ[2][k]*Fc*Gc*aux1;
  }
  //printf("Saiu Hexahedral::mass(%d, %d)\n",m1,m2);
  return aux;
};
//01/07/2014
// ************************************************************************
// Creates the local matrices using static condensation
// It uses NEWMAT matrices and operators
// ************************************************************************
void Hexahedral::make_local_matrices()
{
	//printf("Hexahedral::make_local_matrices\n");
  int ni; 
  int i,ii,j,jj,k;
  ni=nn-nb;
#ifdef PRINTF_ON
  printf("Making local matrices 1 Hexahedral: nb = %d  ni= %d\n", nb,ni);
#endif
  double M[nn][nn];
  int p1,q1,r1,p2,q2,r2;
  // int flag1, flag2;
  double aux, aux0,aux1, Fa,Fb,Fc,Ga,Gb,Gc;
  int m1,m2;
  for(m1=0;m1<nn;m1++){
    // indexes of point m1
    p1= mode_[m1].p_val();
    q1= mode_[m1].q_val();
    r1= mode_[m1].r_val();
    for(m2=m1;m2<nn;m2++){
      // indexes of point m2
      p2= mode_[m2].p_val();
      q2= mode_[m2].q_val();
      r2= mode_[m2].r_val();
      aux=0.0;
      for(k=0;k<Q[2];++k){
	aux0=0.0;
	for(j=0;j<Q[1];++j){
	  aux1=0.0;
	  for(i=0;i<Q[0];++i){
	    Fa=Psia(P[0], p1,xGQ[0][i]);
	    Ga=Psia(P[0], p2,xGQ[0][i]);
	    aux1+=wGQ[0][i]*Fa*Ga;
	  }
	  Fb=Psia(P[1],q1, xGQ[1][j]);
	  Gb=Psia(P[1],q2, xGQ[1][j]);
	  aux0+=wGQ[1][j]*Fb*Gb*aux1;
	}
	Fc=Psia(P[2],r1, xGQ[2][k]);
	Gc=Psia(P[2],r2, xGQ[2][k]);
	aux+=wGQ[2][k]*Fc*Gc*aux0;
      }
      M[m1][m2]=aux;
      M[m2][m1]=aux;
    }
  }
  
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
  printf("Saindo Hexahedral::make_local_matrices: nb = %d  Ni= %d\n", nb,ni);
#endif
  
};

// ************************************************************************

void Hexahedral::make_mass_matrices(int NFields)
{
  //printf("Hexahedral::make_mass_matrices\n");
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
// ************************************************************************
// Estabelece os parametros de Gauss para o quadrilatero
// em coordenadas eta1 e eta2
// (Gauss-Jacobi quadratura) Ver Apendice B de Karniadakis & Sherwin
// ************************************************************************
void Hexahedral::gauss_parameters_default()
{
  //printf("Hexahedral::gauss_parameters\n");
  double x[MAXQ];
  double wtemp[MAXQ];
  double Dtemp[MAXQ][MAXQ];
  
  // Valores para as duas dimensoes
  for(int i=0;i<ndim;i++){ 
    if(gqt[i]==1)
      {//gqt[i]=1 Gauss-Jacobi
	Gauss_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
      }
    
    else if(gqt[i]==2){
      // *******************************\/ \/\/\/\/************************
      Gauss_Radau_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
    }
    
    else{//gqt[i]=3 
      Gauss_Lobatto_Jacobi_parameters(Q[i],0.0,0.0,x,wtemp,Dtemp);
    }
    
    for(int j=0; j<Q[i]; j++){
      xGQ[i][j]=x[j];
      wGQ[i][j]=wtemp[j];
    }
    
    for(int k=0;k<Q[i];k++){
      for(int l=0;l<Q[i];l++){
	D[k][l][i]=Dtemp[k][l];
	//printf(" Hexahedral D[%d][%d][%d] = %g\n", k,l,i,D[k][l][i]);
      }
    }
  }
};
//06/03/2008
// ************************************************************************
//
// ************************************************************************
void Hexahedral::vector_of_integral_of_f_Phi_dv(double vec[],
				      double (*func)(double,double,double),
				      const Vertice vert[],const int map[],
				      const double JV[])
{
  int n,p,q,r;
  int i,j,k;
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  // double e1p,e1m,e2p,e2m,e3p,e3m;
  double x1,x2,x3;//y1,y2,y3,z1,z2,z3;
  int Pdim=P[0]+1;
  double fp[Pdim][Q[1]][Q[2]], fpq[Pdim][Pdim][Q[2]];
  // primeira matriz fp(xi2,xi3)
  for(p=0;p<Pdim;p++){
    // fazer loop sobre xi2 e xi3 
    for(k=0;k<Q[2];k++){
      eta3=xGQ[2][k];
      for(j=0;j<Q[1];j++){
	eta2=xGQ[1][j];
	// somar sobre xi1 (i)
	aux=0.0;
	for(i=0;i<Q[0];i++){
	  // coordenadas
	  eta1=xGQ[0][i];
	  int q = i+Q[0]*(j+Q[1]*k);
	  x1=0;
	  x2=0;
	  x3=0;
	  for(int n=0;n<nv;++n){
	    x1+=vert[map[n]].x * Phi_val[n][q];
	    x2+=vert[map[n]].y * Phi_val[n][q];
	    x3+=vert[map[n]].z * Phi_val[n][q];
	  }
    
	  Fa=Psia(P[0], p,eta1);
	  
	  aux+=(wGQ[0][i]*Fa*func(x1,x2,x3)*JV[q]);// <==== Jacobian JV
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
	  Fb=Psia(P[1],q,eta2);
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
      Fc=Psia(P[2],r,eta3);
      aux+=wGQ[2][k]*Fc*fpq[p][q][k];
    }
    vec[n]=aux;
  }
}
//01/07/2014
// ************************************************************************
//
// ************************************************************************
void Hexahedral::vector_of_integral_of_f_Phi_dv(double vec[],
				      const double func[],
				      const double JV[])
{
  int n,p,q,r;
  int i,j,k;
  double aux,Fa,Fb,Fc,eta1,eta2,eta3;
  //double e1p,e1m,e2p,e2m,e3p,e3m;
 // double x1,x2,x3;//y1,y2,y3,z1,z2,z3;
  int Pdim=P[0]+1;
  double fp[Pdim][Q[1]][Q[2]], fpq[Pdim][Pdim][Q[2]];
  // primeira matriz fp(xi2,xi3)
  for(p=0;p<Pdim;p++){
    // fazer loop sobre xi2 e xi3 
    for(k=0;k<Q[2];k++){
      eta3=xGQ[2][k];
      for(j=0;j<Q[1];j++){
	eta2=xGQ[1][j];

	// somar sobre xi1 (i)
	aux=0.0;
	for(i=0;i<Q[0];i++){
	  // coordenadas
	  eta1=xGQ[0][i];
	  int q = i+Q[0]*(j+Q[1]*k);
	/*  x1=0.0;
	  x2=0.0;
	  x3=0.0;
	  for(int n=0;n<nv;++n){
	    x1+=vert[map[n]].x * Phi_val[n][q];
	    x2+=vert[map[n]].y * Phi_val[n][q];
	    x3+=vert[map[n]].z * Phi_val[n][q];
	  }
    */
	  Fa=Psia(P[0], p,eta1);
	  
	  aux+=(wGQ[0][i]*Fa*func[q]*JV[q]);// <==== Jacobian JV
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
	  Fb=Psia(P[1],q,eta2);
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
      Fc=Psia(P[2],r,eta3);
      aux+=wGQ[2][k]*Fc*fpq[p][q][k];
    }
    vec[n]=aux;
  }
}
//01/07/2014

// ************************************************************************
void Hexahedral::printtofile(FILE * fout,const double u[],
			     double (*func)(double,double,double), 
			     const Vertice vert[], const int map[])
{
  //printf("Hexahedral::printtofile\n");
  double aux,x1,x2,x3;
  double ftemp[Q[0]*Q[1]*Q[2]];
  int m1=0;
  evalGQ(ftemp,u);
  for(int k=0;k<Q[2];++k){
    for(int j=0;j<Q[1];j++){
      for(int i=0;i<Q[0];i++){
        int q = i+Q[0]*(j+Q[1]*k);
        x1=0.0;
        x2=0.0;
        x3=0.0;
        for(int n=0;n<nv;++n){
          x1+=vert[map[n]].x * Phi_val[n][q];
          x2+=vert[map[n]].y * Phi_val[n][q];
          x3+=vert[map[n]].z * Phi_val[n][q];
        }
        aux=ftemp[m1++];
        fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
      }
    }
  }
};
// ************************************************************************
void Hexahedral::printtofile(FILE * fout,const double u[],
				const Vertice vert[], const int map[])
{
	//printf("Hexahedral::printtofile1\n");
  double aux,x1,x2,x3;
  double ftemp[Q[0]*Q[1]*Q[2]];
 
  int m1=0;
  evalGQ(ftemp,u);
  for(int k=0;k<Q[2];++k){
    for(int j=0;j<Q[1];j++){
      for(int i=0;i<Q[0];i++){
    	int q = i+Q[0]*(j+Q[1]*k);
        x1=0.0;
        x2=0.0;
        x3=0.0;
        for(int n=0;n<nv;++n){
          x1+=vert[map[n]].x * Phi_val[n][q];
          x2+=vert[map[n]].y * Phi_val[n][q];
          x3+=vert[map[n]].z * Phi_val[n][q];
        }
        aux=ftemp[m1++];
        fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
      }
    }
  }
};
// ************************************************************************
void Hexahedral::printGQtofile(FILE * fout,const double ftemp[],
				  const double ftemp1[],
				  const Vertice vert[], const int map[])
{
	//printf("Hexahedral::printGQtofile\n");
  double aux, aux1,x1,x2,x3;
  int m1=0;
  for(int k=0;k<Q[2];++k){
    for(int j=0;j<Q[1];j++){
      for(int i=0;i<Q[0];i++){
     	int q = i+Q[0]*(j+Q[1]*k);
        x1=0.0;
        x2=0.0;
        x3=0.0;
        for(int n=0;n<nv;++n){
          x1+=vert[map[n]].x * Phi_val[n][q];
          x2+=vert[map[n]].y * Phi_val[n][q];
          x3+=vert[map[n]].z * Phi_val[n][q];
        }
        aux=ftemp[m1];
        aux1=ftemp1[m1++];
        fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,aux1);
      }
    }
  }
};
// ************************************************************************
// ************************************************************************
void Hexahedral::printwGQtofile(FILE * fout,
				   const Vertice vert[],
				   const int map[],
				   const double JV[])
{
	//printf("Hexahedral::printwGQtofile\n");
  double aux,x1,x2,x3;
  int m1=0;
  for(int k=0;k<Q[2];++k){
    for(int j=0;j<Q[1];j++){
      for(int i=0;i<Q[0];i++){
      	int q = i+Q[0]*(j+Q[1]*k);
	x1=0.0;
	x2=0.0;
	x3=0.0;
	for(int n=0;n<nv;++n){
	  x1+=vert[map[n]].x * Phi_val[n][q];
	  x2+=vert[map[n]].y * Phi_val[n][q];
	  x3+=vert[map[n]].z * Phi_val[n][q];
	}
	aux=wGQ[0][i]*wGQ[1][j]*wGQ[2][k]*JV[m1++];
	fprintf(fout,"%11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux);
      }
    }
  }
};
// ************************************************************************
void Hexahedral::printtoarray(const double u[], 
			    const Vertice vert[], const int map[],
			    double x[], double y[], double z[], double ftemp[])
{
  //printf("Hexahedral::printtoarray\n");
  int m1=0;
  evalGQ(ftemp,u);
  for(int k=0;k<Q[2];++k){
    for(int j=0;j<Q[1];j++){
      for(int i=0;i<Q[0];i++){
      	int q = i+Q[0]*(j+Q[1]*k);
	x[m1]=0.0;
	y[m1]=0.0;
	z[m1]=0.0;
	for(int n=0;n<nv;++n){
	  x[m1]+=vert[map[n]].x * Phi_val[n][q];
	  y[m1]+=vert[map[n]].y * Phi_val[n][q];
	  z[m1]+=vert[map[n]].z * Phi_val[n][q];
	}
	m1++;
      //fprintf(fout,"%11.4e %11.4e %11.4e %11.4e %11.4e\n",x1,x2,x3,aux,func(x1,x2,x3));
      }
    }
  }
};
// ************************************************************************
// Calcular os valores da funcao nos nos de Gauss (overloaded)
// NF = numero de campos (default = 1)
// nvar= numero da variavel (default=0)
// ************************************************************************
void Hexahedral::evalGQ(double f0[],const double u0[],
			const int NF,const int nvar)
{
  //printf("Hexahedral::evalGQ1\n");
  double aux,eta1,eta2,eta3;
  int i,j,k,p,q,r;
  //int Pdim=P[0]+1;
  double fpq[ P[0]+1 ][ P[1]+1 ][ Q[2] ], fp[ P[0]+1 ][ Q[1] ][ Q[2] ];
  //printf("Hexahedral::evalGQ    NF = %d    nvar= %d\n",NF,nvar);
 
  for (k=0;k<Q[2];++k) {
    eta3=xGQ[2][k];
    for(p=0;p<=P[0];++p) {
      for(q=0;q<=P[1];++q) {

	aux=0.0;
	for(r=0;r<=P[2];++r) {
	  aux += Psia(P[2],r,eta3) * u0 [ ind_mode_[p][q][r] ];
	}
	fpq[p][q][k] = aux;
      }
    }
  }
  
  for (j=0;j<Q[1];++j) {
    eta2=xGQ[1][j];
    for (k=0;k<Q[2];++k) {
      for(p=0;p<=P[0];++p) {
	aux=0.0;
	for(q=0;q<=P[1];++q) {
	  aux += Psia(P[1],q,eta2) * fpq [p][q][k];
	}
	fp[p][j][k] = aux;
      }
    }
  }

  for(i=0; i<Q[0]; ++i){
    eta1=xGQ[0][i];
    for(j=0;j<Q[1];j++){
      for (k=0;k<Q[2];++k){
	aux=0.0;
	for(p=0;p<=P[0];++p){
	  aux += Psia(P[0],p,eta1) * fp[p][j][k];
	}
	f0[ NF*(i+Q[0]*(j + Q[1] * k)) + nvar] = aux;
      }
    }
  }
};
//02/07/2014
// ************************************************************************
// Calcular os valores da funcao nos pontos de Gauss
// ************************************************************************
void Hexahedral::evalGQ(double f0[],double f1[],
			const double uh0[],const double uh1[])
{
  //printf("Hexahedral::evalGQ\n");
  double aux0,aux1, Fa,Fb,eta1,eta2;
  int a;
  int i,j,p,q,k;
  int n=0;
  int Pdim=P[0]+1;	
  double ftemp0[Pdim],ftemp1[Pdim];
  for(j=0;j<Q[1];j++){
    eta2=xGQ[1][j];
    // ********************************************************************
    // Construcao dos vetores temporarios
    // ********************************************************************
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
    // ********************************************************************
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
// Evaluates the value of the field at the vertices of the element
// ************************************************************************
void Hexahedral::computeVertice(double f_vert[],const double u[], 
		    const Vertice vert[], const int map[])
{
  //printf("Hexahedral::computeVertice\n");
  double aux, Fa,Fb,eta1,eta2;//x1,x2,x3;
  // double eaux,faux,gaux,haux;
  // double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
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
void Hexahedral::computeAtPoints(const int npoints,
				 const double LocCoord[],
				 const double u[],
				 const Vertice vert[],
				 const int map[],
				 double f[],
				 double GloCoord[])
{
	//printf("Hexahedral::computeAtPoints\n");
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
void Hexahedral::computeFuncGQ(double f_[], 
			       const Vertice vert[], const int map[],
			       double (*func)(double,double,double))
{
  double x[3][8];
  double eta[3];
  double f[2][3]; // funcoes interpolantes
  double aux,aux0,aux1,aux2;
  int a[3][8] = {
    {0,1,0,1,0,1,0,1},
    {0,0,1,1,0,0,1,1},
    {0,0,0,0,1,1,1,1}
  };
  
  for(int i=0;i<8;++i)
  {
    x[0][i]=vert[map[i]].x;
    x[1][i]=vert[map[i]].y;
    x[2][i]=vert[map[i]].z;
  }
  int n=0;
  for(int k=0;k<Q[2];++k){
    eta[2]=xGQ[2][k];
    f[0][2] = 0.5*(1-eta[2]); // funcao atras
    f[1][2] = 0.5*(1+eta[2]); // funcao a frente
    for(int j=0;j<Q[1];++j){
      eta[1]=xGQ[1][j];
      f[0][1] = 0.5*(1-eta[1]); // funcao inferior
      f[1][1] = 0.5*(1+eta[1]); // funcao superior
      for(int i=0;i<Q[0];++i){
        eta[0]=xGQ[0][i];
        f[0][0] = 0.5*(1-eta[0]); // funcao esquerda
        f[1][0] = 0.5*(1+eta[0]); // funcao direita
  
        aux0=0.0;
        aux1=0.0;
        aux2=0.0;
        for(int l=0;l<8;++l){
          aux = f[a[0][l]][0] * f[a[1][l]][1] * f[a[2][l]][2];
          aux0 += x[0][l] * aux;
          aux1 += x[1][l] * aux;
          aux2 += x[2][l] * aux;
        }
        f_[n++]=func(aux0,aux1,aux2);
      }
    }
  }
};
// ***************************************************************************/
//06/03/2008
/*
void Hexahedral::make_Phi(const int m,double Phi[])
{
	//printf("Hexahedral::make_Phi\n");
  //printf("Hexahedral::make_Phi, Q[0] = %d  Q[1]= %d, P[0] = %d, P[1] = %d\n",Q[0],Q[1],P[0],P[1]);
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
void Hexahedral::eval_Phi(const int n,double Phi[])
{
  //printf("Hexahedral::eval_Phi\n");
  int q2=Q[0]*Q[1]*Q[2];
  for(int i=0;i<q2;++i)Phi[i]=Phi_val[n][i];
};
// ************************************************************************
void Hexahedral::eval_GradPhi(const Vertice vert[], const int map[],const int n,double ** der)
{
	//printf("Hexahedral::eval_GradPhi\n");
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


// ************************************************************************
void Hexahedral::Jacobian(const Vertice vert[],const int map[],double JV[])
{
  //cout << "Hexahedral::Jacobian" << endl;
  //int q2=Q[0]*Q[1];
  double x[3][8];
  double eta[3];
  double f[2][3]; // funcoes interpolantes
  double s[2] = { -0.5, 0.5}; // derivadas
  NEWMAT::Matrix A(3,3);
  double aux,aux1;
  int a[3];
  int n;
  
  for(int i=0;i<8;++i)
  {
    x[0][i]=vert[map[i]].x;
    x[1][i]=vert[map[i]].y;
    x[2][i]=vert[map[i]].z;
  }
  
  n=0;
  for(int k=0;k<Q[2];++k){
    eta[2]=xGQ[2][k];
    f[0][2] = 0.5*(1-eta[2]);
    f[1][2] = 0.5*(1+eta[2]);
    for(int j=0;j<Q[1];++j){
      eta[1]=xGQ[1][j];
      f[0][1] = 0.5*(1-eta[1]);
      f[1][1] = 0.5*(1+eta[1]);
      for(int i=0;i<Q[0];++i){
      eta[0]=xGQ[0][i];
        f[0][0] = 0.5*(1-eta[0]);
        f[1][0] = 0.5*(1+eta[0]);
       
        for(int m=0;m<3;++m) {
          for(int n=0;n<3;++n) {
            aux=0.0;
            for(int l=0;l<8;++l){// loop sobre os vertices
              a[0] = (l&1);
              a[1] = (l&2)/2;
              a[2] = (l&4)/4;
              aux1=1.0;
              for(int r=0;r<3;++r){
                if(r==n) aux1 *= s[a[r]];
                else aux1 *= f[a[r]][r];
              }
              aux += x[m][l]*aux1;
            }
           // printf(" elemento A(%d ,%d) = %g\n",m,n,aux);
            A.element(m,n)=aux;
          }
        }
        aux=A.Determinant();
        //printf("Jacobiano = %g\n",aux);
        if(aux<=0.0)printf("Hexahedral::Jacobian Erro: Jacobiano nao-positivo = %g\n",aux);
        JV[n++]=aux;
      }
    }
  }
};
// 01/09/2014

// ************************************************************************
// Calculates the gradient wrt x1, x2
// grad[i][m] = component i of the gradient at the point m(pq)
// ************************************************************************

void Hexahedral::Gradiente(double * grad[],
                              const  double fvec[],
                              const Vertice vert[], const int map[])
{
	//printf("Hexahedral::Gradiente 1\n");
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

void Hexahedral::Gradiente(FILE * fout, double * grad[],
			      const  double fvec[], 
			      const Vertice vert[], const int map[])
{
	//printf("Hexahedral::Gradiente 2\n");
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,eta1,eta2;
  double x3=0.0;
  double a2p,a2m,a1p,a1m;
  double x1,x2;//x12,y1,y2,y12;
  int m;
 // double df[MAXQ*MAXQ][3];
  //double a11,a12,a21,a22, J2D;
  //double b[2][2];
 // double aux0,aux1;
  
  // printf("Calculo do Gradiente de um vetor. Elemento Hexahedral\n");
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

void Hexahedral::Gradiente(FILE * fout, double * grad[],
			      double (*func)(double, double, double), 
			      const Vertice vert[], const int map[])
{
	//printf("Hexahedral::Gradiente 3\n");
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
  //printf("Calculo do Gradiente de uma funcao. Elemento Hexahedral\n");
  Gradiente(fout,grad,fvec,vert,map);
};

void Hexahedral::print_nome(FILE * fout)
{
	//printf("Hexahedral::print_nome\n");
  fprintf(fout,"ELEMENTO QUADRILATERAL\n");
}
// ****************************************************************************
// 20/03/2008

// ************************************************************************
// ************************************************************************
void Hexahedral::Dirichlet(const int face_num,
                            const Vertice vert[],
                            const int vert_map[], // numero global dos vertices dos nos
                            const int nmap[], // mapeamento global dos modos locais
                            const int sgn[],
                            int bflag[],
                            double Xbc[],
                            double (*func)(double,double,double))
{
  //cout << "Entrou em Hexahedral::Dirichlet para a face "<< face_num << "\n";
  const int varn = 0;
  //const int NFields = 1;
  int lvert_map[4]={0,1,2,3};
  int quadmap[MAXMODES];// mapeia modos do quadrilatero sobre os do hexaedro
  int p0,q;
  //int p1;
  int i,ii;
  // construcao do mapa dos vertices locais em vertices globais
  for(int i=0;i<4;++i) {
    lvert_map[i]=vert_map[face[face_num][i]];
    //cout << "lvert_map["<< i << "] = "<< lvert_map[i] << endl;
  }
  q =Q[fd0[face_num]];
  p0=P[fd0[face_num]];
 // p1=P[fd1[face_num]];
 
  //dados do quadrilatero (stdel)
  Quadrilateral * quad = new Quadrilateral(p0,q);
  PhElem<1> * localphel = new PhElem<1>();
  localphel->set_ptr_stdel(quad,quad);
  localphel->set_type(3);// quadrilatero
  localphel->set_ptvert(vert); // array de vertices recebido nos argumentos
  localphel->set_Vert_map(4,lvert_map);
  //localphel->compute_JV(0);// Jacobiano J deve ser calculado antes do vetor
  localphel->inicia_vetores();// zera o vetor b de localphel
  int count = 0;
  localphel->inicia_gbnmap(0,count);// inicia gbnmap[i]=i,sgn[i]=1, i=0,nn-1;
  // Encontrar os nos.
  //Criar quadmap[]; mapeamento dos nos do quadrialtero standard sobre
  // os nos globais
  localFaceModeMap(face_num,quadmap); // mapeia o quadrilatero local no hexaedro
  //cout << "ANTES DE COMPUTE_JV" << endl;
 
  //cout << "Dentro do Hexahedral::Dirichlet" << endl;
  localphel->projetar_C0(nullptr, func, varn);
  //cout << "Dentro do Hexahedral::Dirichlet terminou localphel->projetar_C0 da face " <<face_num<< endl;
  //mapeia os resultados sobre os modos globais
  for(i=0;i<quad->nn_val();i++){
    ii=quadmap[i]; // indice do modo no tetraedro
    Xbc[nmap[ii]]=localphel->show_u(varn,i)*sgn[ii];
    bflag[nmap[ii]]=0; // conhecido
  }
  
  localphel->finaliza_vetores();
  delete localphel;
  delete quad;
 
  //cout << "Saindo Hexahedral::Dirichlet para a face "<< face_num <<endl<<endl;
};

/*
// ****************************************************************************
void Hexahedral::Dirichlet(const int aresta,
                           const Vertice vert[],
                           const int vert_map[],
                           const int gbnmap[],
                           const int sgn[],
                           int bflag[],
                           double X[],
                           double (*f)(double,double,double))
{
  //printf("Hexahedral::Dirichlet\n");
  // cout << "Entrou em Hexahedral::Dirichlet\n";
  int flag=0;
  // **************************************************************************
  // flag = 0 : Dirichlet, valor conhecido, bflag=0
  //      = 1 : valor desconhecido, bflag=1
  // **************************************************************************
  int i,j,k,p,q,nd,ii;
  double xa,ya,xb,yb,za,zb;
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
    temp=gbnmap[0];
    X[temp]=sgn[0]*a0;
    bflag[temp]=flag;
  //temp=gbnmap[1]*NFields+varn;
    temp=gbnmap[1];
    X[temp]=sgn[1]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++){
      ii=i+3;
      //temp=gbnmap[ii]*NFields+varn;
      temp=gbnmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  else if(aresta==1){
    //temp=gbnmap[1]*NFields+varn;
    temp=gbnmap[1];
    X[temp]=sgn[1]*a0;
    bflag[temp]=flag;
    //temp=gbnmap[2]*NFields+varn;
    temp=gbnmap[2];
    X[temp]=sgn[2]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++) {
      ii=i+P[0]+2;
      //temp=gbnmap[ii]*NFields+varn;
      temp=gbnmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  else if(aresta==2){
    //temp=gbnmap[3]*NFields+varn;
    temp=gbnmap[3];
    X[temp]=sgn[3]*a0;
    bflag[temp]=flag;
    //temp=gbnmap[2]*NFields+varn;
    temp=gbnmap[2];
    X[temp]=sgn[2]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++) {
      ii=i+P[0]+P[1]+1;
      //temp=gbnmap[ii]*NFields+varn;
      temp=gbnmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
  else{
    //temp=gbnmap[0]*NFields+varn;
    temp=gbnmap[0];
    X[temp]=sgn[0]*a0;
    bflag[temp]=flag;
    //temp=gbnmap[3]*NFields+varn;
    temp=gbnmap[3];
    X[temp]=sgn[3]*ap;
    bflag[temp]=flag;
    for(i=1;i<p;i++) {
      ii=i+2*P[0]+P[1];
      //temp=gbnmap[ii]*NFields+varn;
      temp=gbnmap[ii];
      X[temp]=sgn[ii]*B.element(i-1);
      bflag[temp]=flag;
    }
  }
};
*/
// 20/03/2008

void Hexahedral::teste(int & v)
{
	//printf("Hexahedral::teste\n");
   v=100000;
}

//Verificar
/*
// ************************************************************************
// Mass Matrix entry M(m1, m2)
// ************************************************************************
// Evaluates the inner product of Phi_m1 by Phi_m2
// ************************************************************************
double Hexahedral::mass(int m1,int m2)
{
 //printf("Hexahedral::mass\n");
  int p,q,r,s;
  int flag1, flag2;
  int j, i;
  double aux, aaux, Fa,Fb,Ga,Gb;
  // printf("\nEntrou Hexahedral::mass(%d, %d)\n",m1,m2);
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
  //printf("Saiu Hexahedral::mass(%d, %d)\n",m1,m2);
 return (aux*JV[0]);// <----------multiplica pelo Jacobiano----<
}; */
//06/03/2008
// *****************************************************************************
// Calcula os tracos de Phi, GradPhi e Jb nos pontos de Gauss sobre as arestas
// ordenando-os de forma que haja coincidencia de pontos dos elementos vizinhos
// *****************************************************************************
void Hexahedral::elem_traces(const Vertice vert[],const int map[],const int sinal[],
			     double *** TP,double **** TGP,
			     double * Jb)
{
  //printf("Hexahedral::elem_traces\n");
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
// ******************************************************************************
void Hexahedral::trace_Jb(const Vertice vert[],const int map[],const int sinal[],
                          double * Jb)
{
    //printf("Hexahedral::elem_traces\n");
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
   // za=vert[map[0]].z;
   // zb=vert[map[1]].z;
   // zc=vert[map[2]].z;
   // zd=vert[map[3]].z;
    
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
void Hexahedral::trace(const int lado,const int qmax,const int sinal,
			  const double *valores,double *saida)
{
  //printf("Hexahedral::trace\n");
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
  
  for(int i=0;i<q;++i){
    saida[i]=valores[ind];
    ind+=inc;
  }
  
  //if(qmax!=q){ //01/10/2013
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
  // ********************************************************
  // Os pontos de Gauss onde estao os tracos sao do tipo
  // Gauss-Jacobi, portanto diferentes dos pontos de Gauss
  // do elemento. Logo toda aresta necessita calcular o traco
  // ********************************************************
  double Jac[qmax],wtemp[qmax],Dtemp[MAXQ][MAXQ];
  Gauss_Jacobi_parameters(qmax,0.0,0.0,Jac,wtemp,Dtemp);
  
  for(int i=0;i<qmax;i++){
    double sum=0.0;
    double y=Jac[i] * sinal;
    for(int k=0;k<q;k++){
      double prod=1.0;
      for(int j=0;j<q;j++){
	if(j!=k)prod*=(y-Old[j]);
      }
      sum+=(temp[k]*prod);
    }
    saida[i]=sum;
  } 
  // 01/10/2013
  //}
  /*
    if(sinal == -1) {// Inverter a ordem dos valores
    for ( int i = 0; i < qmax ; i++ ) temp[i] = saida [ (qmax-1-i) ];
    for ( int i = 0; i < qmax ; i++ ) saida[i] = temp[i];
    }  
  */
};
// revisado em 25/10/2011
const int Hexahedral::aresta_lvert(const int & i, const int & j) const {return aresta[i][j];};
const int Hexahedral::face_lvert(const int & i, const int & j) const {return face[i][j];};
const int Hexahedral::show_nvf(const int &i) const {return nvf[i];};
const int Hexahedral::show_face_tipo(const int &i) const {return face_tipo[i];};
const int Hexahedral::show_fd0(const int &i) const {return fd0[i];};
const int Hexahedral::show_fd1(const int &i) const {return fd1[i];};
const int Hexahedral::show_fv2(const int &i) const {return fv2[i];};
const int Hexahedral::show_ind_mode(const int & i, const int & j, const int & k ) const {return ind_mode_[i][j][k];};
void Hexahedral::superficie_externa(const int Vert_map[],const Vertice vert[],
                                    const int & num_local,
                                    double & area,double normal[3])
{
  int v0,v1,v2;
  v0=Vert_map[face[num_local][0]];
  v1=Vert_map[face[num_local][1]];
  v2=Vert_map[face[num_local][3]]; // Atencao: ultimo indice da sequencia
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
// **********************************************************************
void Hexahedral::localFaceModeMap(const int fnum, int quadmap[])
{
  int dir[3] = {fd0[fnum],fd1[fnum],3-fd0[fnum]-fd1[fnum]};
  int ind[3];
  ind[dir[2]]= fv2[fnum] == 1 ? P[dir[2]] : 0;
  int a = 0;
  
  ind[dir[0]]=0;
  ind[dir[1]]=0;
  quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];// A
  
  ind[dir[0]]=P[dir[0]];
  ind[dir[1]]=0;
  quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];// B
  
  ind[dir[0]]=P[dir[0]];
  ind[dir[1]]=P[dir[1]];
  quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];// C
  
  ind[dir[0]]=0;
  ind[dir[1]]=P[dir[1]];
  quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];// D

  
  // AB
  for(int i=1;i<P[dir[0]];++i){
    ind[dir[0]]=i;
    ind[dir[1]]=0;
    quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];
  }
  
  // BC
 
  for(int j=1;j<P[dir[1]];++j){
    ind[dir[0]]=P[0];
    ind[dir[1]]=j;
    quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];
  }
  
  // DC
  
  for(int i=1;i<P[dir[0]];++i){
    ind[dir[0]]=i;
    ind[dir[1]]=P[dir[1]];
    quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];
  }
  
  // AD
 
  for(int j=1;j<P[dir[1]];++j){
    ind[dir[0]]=0;
    ind[dir[1]]=j;
    quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];
  }
  
  //Verificar
  // Interior modes : j runs fastest
  for(int i=1; i<P[dir[0]]; ++i){
    ind[dir[0]]=i;
    for(int j=1; j<P[dir[1]];++j){
      ind[dir[1]]=j;
      quadmap[a++]=ind_mode_[ind[0]][ind[1]][ind[2]];
    }
  }
  
#ifdef PRINTF_ON
  printf("quadmap da face %d\n",fnum);
  for(int i=0;i<a;++i){
    printf("quadmap[%d] = %d ", i,quadmap[i]);
    mode_[quadmap[i]].print();
  }
  printf("Saindo de Hexahedral::localFaceModeMap\n");
#endif
  
};
