
# include "spectral.h"
void quad_ordem(const int Av[4],const int Bi[4],int dir[3], int sgn[3], int & v2)
{
    // Entrada Av na ordem anti-horaria
    //    2
    //    D ----- C 3   Av = {A, B, C, D} nos globais
    //    |       |     Bi = {0, 1, 3, 2} numeracao local da face 0 do Hexahedro usado como exemplo
    //    |       |     Observe que a numeracao local do Hexahedro e peculiar de modo a determinar
    //    A ----- B     a direcao da aresta pela diferenca dos numeros locais dos vertices
    //    0       1
    // Bi  contem os dados da face a ser padronizada.
    //
  // Coloca os vertices do quadrilatero na sequencia padrao
  // ordena os indices de acordo com os numeros dos vertices
  // Encontra o valor minimo e seu indice
  int v[4];
  int i[4];
  int i_min = 0;  // ponto inicial
  int v_min = Av[0]; // valor inicial
  for(int k=1;k<4;++k) {
    if( Av[k] < v_min ) {
      v_min = Av[k];
      i_min = k;
    }
  }
  //Girar
  int count = 0;
  for(int k=0;k<4;++k) {
    count = (k + i_min) % 4;
    v[k] = Av[count];
    i[k] = Bi[count];
  }
  int temp;
  if(v[1] > v[3] ) {
    temp = i[1];
    i[1] = i[3];
    i[3] = temp;
    temp = v[1];
    v[1] = v[3];
    v[3] = temp;
  }
  printf("\n");
  for(int k=0; k < 4; ++k)
    printf("v[%d] = %d  i[%d] = %d\n", k, v[k], k, i[k]);

  // Usar este trecho no hexahedral.cc
  dir[0]=(i[1]-i[0]);

  dir[1]=(i[3]-i[0]);

  if(dir[0] < 0) {sgn[0] = -1;  dir[0] *= sgn[0];} else sgn[0] = 1;
  if(dir[1] < 0) {sgn[1] = -1;  dir[1] *= sgn[1];} else sgn[1] = 1;


  dir[2]=((dir[0] ^ 7) ^ dir[1]);

  v2=(i[0] & i[1] & i[2] & i[3] & dir[2])/dir[2];

  dir[0] /=2;
  dir[1] /=2;
  dir[2] /=2;

  printf("quad_ordem\ndir[0] = %d sgn[0] = %d\ndir[1] = %d sgn[1] = %d\ndir[2] = %d v2 = %d\n",dir[0], sgn[0], dir[1],sgn[1],dir[2], v2);
};


// *********************************************************************
// ordena seq de n inteiros de acordo com a base A
// reordena A e seq simultaneamente de modo que A fique crescente
// retorna em A a base ordenada e em seq a sequencia equivalente
// *********************************************************************
void ordenar(const int n, int A[],int seq[])
{
	int temp, itemp;
	for(int i=0; i<n-1;++i) {
		temp=A[i];
		itemp=seq[i];
		for(int j=i+1;j<n;++j) {
			if(A[j]<temp){
				A[i]=A[j];
				A[j]=temp;
				temp=A[i];
				seq[i]=seq[j];
				seq[j]=itemp;
				itemp=seq[i];
			}
		}
	}
};

// ***********************************************************
double norma_max(int N,double V[])
{
  double norma=0.0;
  double aux,aux1;
  for(int i=0;i<N;i++){
    aux1=V[i];
    aux = (aux1>0) ? aux1 : -aux1;
    if( aux > norma) norma=aux;
  }
  return norma;
}
// ***********************************************************
double norma_l2(int N,double B[])
{
  double valor0=0.0;
  for(int iN=0;iN<N;iN++) valor0 += B[iN]*B[iN];
  valor0 = sqrt(valor0);
  return valor0;
}

int Sinal(int i) // para elemento linear
{
  if(i==0) return (1);
  else return (-1);
};
// ***********************************************************
double funcao(double x, double y,double z)
{
 // double xx = x;
  //double yy = y;
  // double zz = z;
  // double aux,value;
  /*
    double PI=atan(1.0)*4.0;
    aux=x+y;
    if(aux<2.2 && aux>2.0) {
    value=0.5*(1.0+sin((aux-2.1)*5.0*PI));
    }
    else if(aux<=2.0)value=0.0;
    else value=1.0;
    return value;
  */
  //x/=100.0;
  /*  aux=(1.0 + x*x*(1.6*x-2.4));
////aux=1.0-x*x*(3.0-2.0*x);
if (x<0.0) aux=1.0;
if(x>1.0) aux=0.2;
return (aux);
  */
  //return(sin(x));
  return (x*x*x + y*y + z);
};
// ***********************************************************
// Derivada
double g(double x, double y,double z)
{
 // double xx = x;
 // double yy = y;
  //double zz = z;
  return (3.0*y*y);
  //return (cos(x));
  /*
    double aux;
    if(x>1.0)aux=0.0;
    else aux=(6.0*x*(x-1.0));
    return (aux);
  */
};


// ****************************************************************************
// Funcao para determinacao do vetor
// ****************************************************************************


//   double force0(double x, double y,double z)
//   {
//     return(0.0);//funcao(x,y,z);
//   };
//   double force1(double x, double y,double z)
//   {
//     return(1.0);//f(x,y,z);
//   };
//
//   double funcao_contorno(double x, double y,double z)
//   {
//     return funcao(x,y,z);
//   };
//

//   // ****************************************************************************
//   double p1(double x, double y,double z)
//   {
//     return (3.0e6);
//   };
//   double p2(double x, double y,double z)
//   {
//     return (1.0e6);
//   };
//   double fluxo_entrando(double x, double y,double z)
//   {
//     return (1.0e-2);
//   };
//void escreva(Vector B)
//{
//  printf("Chamada de funcao passando Vector como argumento\n%11.4e\n",B(0));
//};

//    // ****************************************************************************

// Testado em 19/07/2014
// *************************************************************************************
// Faz a mascara de faces de antiga para nova (face_mask) para o tetrahedro re-ordenado
// face_mask[num da face original] = num nova face (no tetraedro ordenado padrao)
// *************************************************************************************
void tetrahedro_faz_face_mask(const int & n_in,int ver_temp[],std::vector<int> & face_mask)
{
  int seq_orig[4] = {0,1,2,3}; // sequencia original eh 0,1,2,3
  ordenar(4,ver_temp,seq_orig); // seq_orig retornada representa a nova ordem de acodo com vert_temp
  int facemap[4][3] = // sequencia das faces do tetraedro. Certifique-se que e igual ao de Tetrahedral.
  {          //  sequencia de faces adotada no Gambit. Copiada de Tetrahedral.h
    {0,1,2}, // ABC
    {0,1,3}, // ABD
    {0,2,3}, // ACD
    {1,2,3}  // BCD
  };
  
  int temp[3],temp1[3] = {0,1,2};
  // face nova i em termos de seq_orig
  for(int i=0; i<4; ++i) { // nova face
    for(int j=0; j<3; ++j) temp[j] =seq_orig[facemap[i][j]]; //nova face = i
    ordenar(3,temp,temp1);
    // comparar se temp se iguala a face m
    int flag=0;
    for(int m=0; m<4 && flag == 0; ++m){ //face antiga = m
      flag=1;
      for(int k=0; k<3 && flag ==1; ++k){
        if(temp[k] != facemap[m][k]) flag = 0; // faces diferem
      }
      if(flag == 1) face_mask[m]=i; //faces coincidem
    }
  }
};

// Hexahedral internal Face Modal Connectivity
// To locate uniquely the modes on the faces of neighboring elements
void Hexahedral_face_modal_connectivity(int Av[4], int Bi[4], int P[],int mode_mask[],int sgn_mask[])
{
  int ini[3], fim[3],sinal[3];
  int inc[3];
  int dir[3];
  int sgn[3]={1,1,1};
  int v2;

  quad_ordem(Av,Bi,dir,sgn,v2);

  int iaux = dir[0];
  ini[iaux]=0;
  fim[iaux]=P[iaux] - 1;
  inc[iaux] = P[dir[1]] - 1;
  sinal[iaux] = sgn[0];

  iaux = dir[1];
  ini[iaux]=0;
  fim[iaux]=P[iaux] - 1;
  inc[iaux] = 1;
  sinal[iaux] = sgn[1];

  inc[dir[2]] = 0;
  sinal[dir[2]] = 1;
  if(v2 == 0) {
    ini[dir[2]] = 0; fim[dir[2]] = 1;
  }
  else if(v2 == 1) {
    ini[dir[2]] = P[dir[2]]; fim[dir[2]] = P[dir[2]] + 1;
  }

  for(int k=0;k<3;++k)
    printf("ini[%d] = %d fim[%d] = %d inc[%d] = %d\n",k,ini[k],k,fim[k],k,inc[k]);
  int count=0;
  int n0 = 10;
  printf("Dentro de hex_face_modal_connect\n");
  int aux0,aux1,aux2;
  aux0=1;
  for(int p=ini[0];p < fim[0];++p){
    aux1=1;
    for(int q=ini[1];q < fim[1];++q){
      aux2=1;
      for(int r=ini[2];r < fim[2];++r){
        mode_mask[count] = p*inc[0] + q*inc[1] + r*inc[2] + n0;
        sgn_mask[count] = aux0 * aux1 * aux2;
        printf("count= %d %d %d\n",count, mode_mask[count], sgn_mask[count]);
        aux2 *= sinal[2];
        ++count;
      }
      aux1 *= sinal[1];
    }
    aux0 *= sinal[0];
  }
  //printf("%d %d\n",count,mode_mask[count]);
};
// ********************************************************************************

// ****************************************************************************
// ****************************************************************************
int get_linha(FILE * finput,char linha[256])
// ****************************************************************************
{
  int i=0;
  char c;
//  static int m=0;
  while((c=getc(finput)) == '\n' || c == '\0');
  /* printf("while: %c", c);*/

  if(c==EOF) {
    return i;
  }
  else {
    linha[0]=c;
    i=1;
    while ((c=getc(finput))!='\n' &&  i < 255 && c != EOF)
      {
        linha[i]=c;
        i++;
      }
    linha[i]='\0';
  }
  //printf("Getlinha %d: %s\n", ++m,linha);
  return i;
};
// Procurar esquema mais simples para numerar os  modos das bordas !!!!
// ****************************************************************************
// ****************************************************************************
int teste_aresta(int p,int n0,int n1,int& sign,
		  int& ng, // candidato a sair
		  int& NG,int& NL, std::vector<EDGE>& border,
		  int Ng[],// candidato a sair
		  int nel,int naresta,
		  const Vertice * vert,const int sinal_normal)
// ****************************************************************************
// ****************************************************************************
{
  int flag,edge_num;
  sign=1;
  //printf("NL = %d\n",NL);
  // Verificar se a aresta ja foi numerada
  flag=0; // reset flag; flag = 0 equivale a aresta nao numerada

  for (int i=0; i<NL && flag==0; i++){
    //printf("NL = %d i=%d\n",NL,i);
    if(border[i].Na==n0 && border[i].Nb==n1){
      //ng= Ng[i];//candidato a sair
      flag=1;
      border[i].num_elem += 1;
      border[i].elemento[1]=nel;
      border[i].num_local[1]=naresta;
      border[i].sinal[1]=1;
      border[i].tipo=2;// aresta interior
      edge_num = i;
    }
    else if(border[i].Na==n1 && border[i].Nb==n0){// aresta ja percorrida em sentido oposto
      //ng= Ng[i];//candidato a sair
      flag=1;
      sign=-1;
      border[i].num_elem += 1;
      border[i].elemento[1]=nel;
      border[i].num_local[1]=naresta;
      border[i].sinal[1]=-1;
      border[i].tipo=2;// aresta interior
      edge_num = i;
    }
  }
  // flag = 1 se a aresta ja foi numerada antes

  if(flag==0){// aresta nova
    EDGE temp_border;
    //ng=NG;//candidato a sair
    // calculo da normal
    temp_border.num_elem = 1;
    temp_border.Na=n0;
    temp_border.Nb=n1;
    temp_border.elemento[0]=nel;
    temp_border.num_local[0]=naresta;
    temp_border.sinal[0]=1;
    temp_border.elemento[1]=0;
    temp_border.num_local[1]=0;
    temp_border.sinal[1]=123;
    temp_border.tipo=0;// aresta de contorno (inicialmente todas sao no-flow)

		if(sinal_normal!=0) { // gira vetor (x,y) de sinal_normal*pi/2 no sentido anti-horario
			double x = vert[n1].x-vert[n0].x;
			double y = vert[n1].y-vert[n0].y;
			double modulo = sqrt(x*x+y*y);
			temp_border.comprimento=modulo;
			temp_border.normal[0] = -sinal_normal*y/modulo;
			temp_border.normal[1] =  sinal_normal*x/modulo;
		}
		else { // Caso 3D onde sinal_normal=0
			double x = vert[n1].x-vert[n0].x;
			double y = vert[n1].y-vert[n0].y;
			double z = vert[n1].z-vert[n0].z;
			double modulo = sqrt(x*x+y*y+z*z);
			temp_border.comprimento=modulo;
			temp_border.normal[0] = x/modulo;
			temp_border.normal[1] = y/modulo;
			temp_border.normal[2] = z/modulo;
		}
    border.push_back(temp_border); // insere copia de temp_border no final de border

    //Ng[NL]=NG;//candidato a sair
    edge_num = NL;
    NL++; // acrescenta o numero de arestas armezandas
    // NG+=(p-1);// acrescenta numero de nos armazenados //candidato a sair
  }

  if(NL > MAXNL || NG > MAXNG){
    printf("Ajuste indices\n");
    printf("NL = %d (MAXNL = %d)\nNG = %d (MAXNG = %d)\n",NL,MAXNL,NG,MAXNG);
  }
  return ( edge_num);
};

// *******************************************************************
int aresta_gbnum(const int na, const int nb, int & NL,std::vector<ARESTA> & aresta_vector)
{
  int num=0;
  int flag=0;
  for(int i=0; i<NL && flag==0; ++i) {
    if(aresta_vector[i].Na==na && aresta_vector[i].Nb==nb){
      num=i;
      flag=1;
    }
  }
  if(flag==0){//aresta nova
    ARESTA temp;
    temp.Na=na;
    temp.Nb=nb;
    aresta_vector.push_back(temp);
    num=NL;
    NL++;
  }
  return (num) ;
};
// ***********************************************************************
int face_gbnum(const int nvf, const int v[], int & NF, std::vector<FACE> & face_vector)
{
 // cout << "entrou face_gbnum " << endl;
  int flag = 0;
  int face_n;
  std:vector<int> var;
  for(int i=0;i<nvf;++i) {
    var.push_back(v[i]);
  }
  std::sort(var.begin(),var.end());
  for(int i=0;i<NF && flag==0; ++i) {
    flag=0;// face_vector[i] difere de var
    if(nvf == face_vector[i]._nv){
      flag=1; // face_vector[i] pode coincidir com var
      for(int j=0;j<nvf;++j){
       if( var[j] != face_vector[i]._vertice[j])
         flag=0; // face_vector[i] difere de var
      }
    }
    if(flag == 1) face_n = i;// face_vector[i] coincide com var
  }

  if(flag == 0){ // face nova
    FACE temp_face;
    temp_face._tipo = 0;
    temp_face._nv = nvf;
    for(int j=0;j<nvf;++j)
      temp_face._vertice.push_back(var[j]);
    face_vector.push_back(temp_face);
    face_n = NF;
    NF++;
  }
  return (face_n);
};

// *******************************************************************
int ResolverSistema(const int NumD, const int count,
  		    int * Ti, int * Tj, double * Tx, double * B, double * x)
// *******************************************************************
{
  std::vector<int>     Ap  (NumD+1);
  std::vector<int>     Ai  (count);
  std::vector<double>  Ax  (count);
  std::vector<int>     Map (count);

  for(int i=0; i<NumD; i++)
    x[i]=0.0;

  umfpack_di_triplet_to_col(NumD,NumD,count,Ti,Tj,Tx,&Ap[0],&Ai[0],&Ax[0],&Map[0]);

  double *null = (double *) NULL ;
  void *Symbolic, *Numeric;
  double Mx, Ex, Info [UMFPACK_INFO] ;
  int status;

  (void) umfpack_di_symbolic(NumD,NumD,&Ap[0],&Ai[0],&Ax[0],&Symbolic, null, null);
  (void) umfpack_di_numeric (&Ap[0],&Ai[0],&Ax[0], Symbolic, &Numeric, null, null) ;
  status = umfpack_di_get_determinant (&Mx, &Ex, Numeric, Info) ;
  umfpack_di_free_symbolic (&Symbolic);
  (void) umfpack_di_solve(UMFPACK_A,&Ap[0],&Ai[0],&Ax[0],x,B, Numeric, null, null);

  umfpack_di_free_numeric (&Numeric) ;
  //printf("Terminou resolucao do sistema global\n");

  double det ;
  det = Mx * pow (10.0, Ex) ;
  printf("Terminou resolucao do sistema global\n determinante %g\n",det);

  return (0);
};
// *********************************************************
// Zerar Matriz A e Vetor B
// *********************************************************
void ZerarMatrizVetor(Teuchos::RCP<Epetra_FECrsMatrix> A,
                      Teuchos::RCP<Epetra_FEVector> RHS)
{
  A->PutScalar(0.0);
  RHS->PutScalar(0.0);
}
// ****************************************************************************
// ****************************************************************************
