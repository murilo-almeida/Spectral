#ifndef _GeProb_headers
#define _GeProb_headers

#include "Funcoes_c.h"

// ****************************************************************************
// Classe GeProb (Generic Problem)
// Encapsula os dados para montar e resolver o problema algebrico generico
// Versao com template gerada em 28/10/2016
// ****************************************************************************

// Flag to tell the evaluate routine what objects to fill
enum FillType {F_ONLY, MATRIX_ONLY, ALL};

template <typename ElemType,int N_VAR=1, int N_FIELDS=1>
class GeProb
{
 public:

  GeProb(Epetra_Comm& comm);
  ~GeProb();
  //Epetra_Comm* show_Comm(){return Comm;};
  Teuchos::RCP<Epetra_Comm> show_Comm(){return Comm;};
  Teuchos::RCP<Epetra_Vector> getSolution(){return solution;};
  // Return a reference to the Epetra_Vector with the Jacobian
  Teuchos::RCP<Epetra_CrsMatrix> getJacobian();


  //----------------------------------------------------------------------
  // Funcoes Gerais de definicao da malha

  // Usada para ler arquivo no formato do Gambit
  // void Ler_Header_Gambit();
  // void Ler_Coordenadas_Gambit();

  void Ler__Geometria(FILE *finput,double *coord,int *Ed,int *BC);

  void set_id(const int rank, const int size){myid=rank;comm_size=size;};
 // void set_dim(const int & d) {dim = d;};
  void Marcar_condicoes_contorno(int *, std::vector < std::vector<int> >);
  void Processar_elementos();

  //void set_bflag(int * ptr);
 // void set_novoNum(int * ptr);

  void set_finput(FILE * ptr);
  void set_fout(FILE * ptr);
  void set_fout1(FILE * ptr);
  void set_fout2(FILE * ptr);
  void set_fout3(FILE * ptr);
  int show_NG();
  int show_NUMNP();
  int show_NELEM();
  int show_NL(){return NL;};
  int show_NF(){return NF;};
  int show_NumD(){return NumD;};
  void RenumerarNos();
  void RenumerarNos(int,int&);
  void Particionar_malha(const int * buf);
  void ResolverComTrilinos(const std::string Pack,
                           Epetra_Map  Map,
                           Teuchos::RCP<Epetra_FECrsMatrix>  A,
                           Teuchos::RCP<Epetra_FEVector>  RHS,
                           double_t & norm_delta_X);
  void Construir_bordas();

  void gbnmap_continuous(const int & ivar,int & count);
  void gbnmap_continuous(int & count);
  void mesh_suporte(const int & P, int Temp_V[], int Temp_A[], int Temp_F[],int & count);
  void gbnmap_discontinuous(const int & ivar,int & count);
  void Facet_eco(FILE * fout);
  void MPI_Recebe_Dados_GeProb();
  void Ler_e_Processar_malha(char *arq_geo);
 // void teste_border_facet();


protected:

    FillType flag;
    //Epetra_Map *StandardMap;
    Teuchos::RCP<Epetra_Map> StandardMap;
    //Epetra_Map *OverlapMap;
    //Epetra_Import *Importer;
    //Epetra_Vector *RHS;
    Teuchos::RCP<Epetra_FEVector>  RHS;
    //Teuchos::RCP<Epetra_FECrsGraph> AA;
    Teuchos::RCP<Epetra_FECrsMatrix> A;
    Epetra_Comm *Comm;

  Teuchos::RCP<Epetra_Vector> solution;


  int MyPID;              // Process number
  int NumProc;            // Total number of processes
  int NumElemTypeents;    // Number of elements owned by this process
  int NumGlobalElements;  // Total Number of elements



  int myid = 0;
  int comm_size = 0;

  // ******************************************
  // Dados das funcoes do problema
  int NumD = 0;
  int NumC = 0;
  const int NumVAR = N_VAR; // Para ser usado nas classes derivadas
  const int NumFIELDS = N_FIELDS;// Para ser usado nas classes derivadas
  Field_struct Field[N_FIELDS];

  int FieldOfVar[N_VAR];

  // ******************************************
  // Dados da Grade
  int dim = 0; // nao iniciada
  int NG = 0;
  int NUMNP = 0;
  int NELEM = 0;
  int NDFCD = 0;
  int NL = 0;
  int NF = 0;
  int NBORDER = 0;
  Vertice * V = nullptr;
  std::vector<ARESTA>   Aresta;
  std::vector<FACE>     Face;
  ElemType * el = nullptr; // Aqui eh o ponto chave de usar o <typename ElemType>
  //Elemento * Elem;
  int * novoNum = nullptr;


  //std::vector<EDGE>     border;
  EDGE * border = nullptr;
  std::vector<FACET>    facet;

  // ******************************************
  // dados das condicoes de contorno
  std::vector<int> bflag;
  int DNBC = 0;
  int NBC[20][3];
  int nin = 0;
  int nout = 0;
  std::vector<int> in_borders, out_borders;

  // ******************************************
  // Dados da particao dos vertices da malha
  int NumPart; //!< Numero de particoes
	std::vector<PARTICAO> Particao;

  // ******************************************
  Linear * ptrLinear[N_FIELDS];
  Triangle * ptrTriang[N_FIELDS];
  Quadrilateral * ptrQuadri[N_FIELDS];
  Tetrahedral * ptrTetrahedral[N_FIELDS];
  Hexahedral *  ptrHexahedral[N_FIELDS];

  // ******************************************
  // Para usar Trilinos

  //std::vector<int> NumNz;//<! numero de nao zeros em todas linhas
  //std::vector<int> NNz;  //<! numero de nao zeros nas linhas pertencentes ao processo
  std::string TrilinosSolver; // tipo de solver usado
  std::string AmesosSolverType; // se usar amesos qual o tipo duas opcoes: Amesos_Umfpack e Amesos_Klu
};
// **************** End of class declaration *********************************

/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
GeProb<ElemType,N_VAR,N_FIELDS>::GeProb(Epetra_Comm& comm): Comm(&comm)
{};

template <typename ElemType,int N_VAR,int N_FIELDS>
GeProb<ElemType,N_VAR,N_FIELDS>::~GeProb()
{cout<<"\nGeProb destructor\n";

    if (V!=nullptr){
  delete [] V; V=nullptr;
    }
    if(el!=nullptr){
  delete [] el; el=nullptr;
    }
    if(novoNum != nullptr){
  delete [] novoNum; novoNum = nullptr;
    }
  // Terminou de estabelecer as condicoes de contorno

  if(border != nullptr) {
    // cout << "liberando border"<< endl;
    for(int i=0;i<NBORDER;++i){
      if(border[i].pdir != nullptr) {
        delete [] border[i].pdir;
        border[i].pdir = nullptr;
       // cout << "liberando border.pdir"<< endl;
      }
      if(border[i].sdir != nullptr) {
        delete [] border[i].sdir;
        border[i].sdir = nullptr;
       // cout << "liberando border.sdir"<< endl;
      }
    }
    delete [] border; border = nullptr;
     cout << "liberou border"<< endl;
  }

 for(int i=0;i<N_FIELDS;++i){
   // cout<< "deletando ponteiros dos elementos padroes; passada "<< i << endl;
    if (ptrLinear[i] != nullptr){delete ptrLinear[i]; ptrLinear[i] = nullptr;}
   // cout << "1" << endl;
    if (ptrTriang[i] != nullptr){delete ptrTriang[i]; ptrTriang[i] = nullptr;}
   // cout << "2" << endl;
    if (ptrQuadri[i] != nullptr){delete ptrQuadri[i]; ptrQuadri[i] = nullptr;}
    //cout << "3" << endl;
    if (ptrTetrahedral[i] != nullptr){delete ptrTetrahedral[i]; ptrTetrahedral[i] = nullptr;}
    //cout << "4" << endl;
    if (ptrHexahedral[i] != nullptr){delete ptrHexahedral[i]; ptrHexahedral[i] = nullptr;}
    //cout << "5" << endl;
  }

    cout << "liberou pointer dos elementos espectrais"<< endl;

    //delete StandardMap;

    cout << "Saiu de GeProb destructor"<< endl;
};
// ****************************************************************************
// ****************************************************************************
//template <typename ElemType,int N_VAR,int N_FIELDS>
//void GeProb<ElemType,N_VAR,N_FIELDS>::set_Vertice_array(Vertice * ptvert)
//{
 // V=ptvert;
//};
/*****************************************************************************/
/*****************************************************************************/

//template <typename ElemType,int N_VAR,int N_FIELDS> void GeProb<ElemType,N_VAR,N_FIELDS>::set_PhElem_array(PhElem * ptr)
//{
//  el=ptr;
//};

/*****************************************************************************/
/*****************************************************************************/
// template <typename ElemType,int N_VAR,int N_FIELDS> void GeProb<ElemType,N_VAR,N_FIELDS>::set_EDGE_array(EDGE *ptr)
// {
//   border=ptr;
// };
///***************************************************************************/
///***************************************************************************/
//template <typename ElemType,int N_VAR,int N_FIELDS>
//void GeProb<ElemType,N_VAR,N_FIELDS>::set_bflag(int * ptr)
//{
//  // bflag=ptr;
//};
///***************************************************************************/
///***************************************************************************/

//template <typename ElemType,int N_VAR,int N_FIELDS>
//void GeProb<ElemType,N_VAR,N_FIELDS>::set_novoNum(int * ptr)
//{
//  novoNum=ptr;
//};

/*****************************************************************************/
/*****************************************************************************/
// *****************************************************************************
// *****************************************************************************

/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::set_finput(FILE * ptr)
{
  // finput=ptr;
};
/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::set_fout(FILE * ptr)
{
  //fout=ptr;
};
/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::set_fout1(FILE * ptr)
{
  //fout1=ptr;
};
/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::set_fout2(FILE * ptr)
{
  //fout2=ptr;
};
/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::set_fout3(FILE * ptr)
{
  //fout3=ptr;
};
/*****************************************************************************/
/*****************************************************************************/
// template <typename ElemType,int N_VAR,int N_FIELDS>
// void GeProb<ElemType,N_VAR,N_FIELDS>::set_X0_ptr(double * ptr)
// {
//   X0=ptr;
// };
/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
int GeProb<ElemType,N_VAR,N_FIELDS>::show_NG()
{
  return NG;
};
/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
int GeProb<ElemType,N_VAR,N_FIELDS>::show_NUMNP()
{
return NUMNP;
};
/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
int GeProb<ElemType,N_VAR,N_FIELDS>::show_NELEM()
{
return NELEM;
};

/*****************************************************************************/
/*****************************************************************************/
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::Processar_elementos()
{
  int ng=NUMNP;// Numero global, inicialmente = NUMNP(number of mode points)
  int type_num;

  // inicializa os elementos padrões
  //int n = N_FIELDS;//Field.size();

  for(int i=0; i < N_FIELDS; ++i) {
    //	cout << "DG_Prob::set_orders:  "<< i << '\n';
    // cout << "Aloca memoria dinamica\n ";
    ptrLinear[i] = new Linear(Field[i].P,Field[i].Q);// type =1
    ptrTriang[i] = new Triangle(Field[i].P,Field[i].Q);// type =2
    ptrQuadri[i] = new  Quadrilateral(Field[i].P,Field[i].Q);// type =3
    ptrTetrahedral[i] = new  Tetrahedral(Field[i].P,Field[i].Q);// type =4
    ptrHexahedral[i] =  new  Hexahedral(Field[i].P,Field[i].Q);// type =5
  }
  //cout << "Passou N_FIELDS:  "<< n << '\n';
  Stdel * ptstdel[N_VAR];

  NL=0;// Numero de arestas (lados)
  NF=0;// Numero de faces

  for(int i=0;i<NELEM;++i){
    type_num=el[i].type_val(); // Tipo do elemento

    switch(type_num){

      case 1:
        for(int j=0;j<N_VAR;++j) ptstdel[j]=ptrLinear[FieldOfVar[j]];
      break;

      case 2:
        for(int j=0;j<N_VAR;++j) ptstdel[j]=ptrTriang[FieldOfVar[j]];
      break;

      case 3:
        for(int j=0;j<N_VAR;++j) ptstdel[j]=ptrQuadri[FieldOfVar[j]];
      break;

    case 4:
        for(int j=0;j<N_VAR;++j) ptstdel[j]=ptrTetrahedral[FieldOfVar[j]];
      break;

      case 5:
        for(int j=0;j<N_VAR;++j) ptstdel[j]=ptrHexahedral[FieldOfVar[j]];
        break;

      default:
        printf("Elemento lido nao eh de tipo conhecido\n");
        exit(0);
    }

    el[i].set_ptr_stdel(ptstdel);

    // ***********************************************************

    // Processamento do elemento
    el[i].Processar_dados(NL,Aresta,NF,Face);
  }
  // *******************************************
  // NG=Numero de nos vezes numero de variveis *
  // NG=ng*N_VAR;
  NG=ng;
  // *******************************************

  if(NL>MAXNL){
    printf("ERRO: Numero de arestas NL(= %d) > MAXNL(=%d) !!!\n",NL,MAXNL);
    printf("Aumente MAXNL !!!\n");
    // exit(0);
  }
};


// 26/11/2013

// ****************************************************************************
// ****************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::RenumerarNos()
// ****************************************************************************
// ****************************************************************************
{

  // **************************************************************************
  // Renumeracao dos graus de liberdade globais
  // Renumera primeiro os desconhecidos (bflag=1)
  // **************************************************************************
  //cout << "template <typename ElemType,int N_VAR,int N_FIELDS> GeProb<ElemType,N_VAR,N_FIELDS>::RenumerarNos()\n";
  //cout << "Aloca memoria dinamica\n";
  novoNum= new int [NG];
  int a=0;
  for(int i=0;i<NG;i++){
    if(bflag[i]==1){
      novoNum[i]=a;
      a++;
    }
  }
  // time_t time1;
  // time(&time1);

  // **************************************************************************
  // Renumera as variaveis conhecidas (bflag=0)
  // **************************************************************************
  NumD=a; // Numero de desconhecidos
  for(int i=0;i<NG;i++){
    if(bflag[i]==0){// Conhecidos
      novoNum[i]=a;
      a++;
    }
  }

#ifdef ECHO_ON
  //printf("Numero Total %d Numero de desconhecidos %d\n", NG,NumD);
  fprintf(fout1,"Echo novoNum\n");
  for(int i=0;i<NG;i++)
    fprintf(fout1,"novoNum[%3d]= %3d   bflag[%3d]= %3d\n",i,novoNum[i],i,bflag[i]);
#endif
  printf("PONTO 6: Renumerou os nos. clock acumulado=%u\n",(unsigned)clock());
}
// ****************************************************************************
// ****************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::RenumerarNos(int numf,int & ND)
// ****************************************************************************
// ****************************************************************************
{
  // **************************************************************************
  // Renumeracao dos graus de liberdade globais
  // Renumera primeiro os desconhecidos (bflag=1)
  // **************************************************************************

  // **************************************************************************
  // Renumera as variaveis conhecidas (bflag=0)
  // **************************************************************************

  int ng = NG/N_VAR;
  int ii;
  int a=0;
  for(int i=0;i<ng;i++){
    ii=i*N_VAR+numf;
    if(bflag[ii]==1){//<=====Desconhecido===<
      novoNum[ii]=a;
      a++;
    }
  }
  ND=a; // Numero de desconhecidos
  for(int i=0;i<ng;i++){
    ii=i*N_VAR+numf;
    if(bflag[ii]==0){//<===Conhecido====<
      novoNum[ii]=a;
      a++;
    }
  }

};


// *****************************************************************************
// *****************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::Particionar_malha(const int * buf)
// *****************************************************************************
{
  if(NumPart > 1) {
		int np;
    int count = 0;
    // Designar a particao dos vertices
    for(int i = 0; i < NumPart; ++i) {
      np=buf[count++];// numero de vertices da particao i
      for(int j = 0; j < np; ++j)
				V[buf[count++]].part_num=i;
    }
  }
  else {
    for(int i = 0; i < NUMNP; ++i)
      V[i].part_num=0;
  }
  // ********************************************
  // Designar a particao dos elementos e bordas
  // ********************************************
  for(int i = 0; i < NELEM; ++i) {
    el[i].set_part_num();
  }
  // *****************************************************************
  // Designar a particao 0 para os elementos com condicoes de contorno
  // de entrada (inflow) ou saida (outflow)
  // *****************************************************************


  // ************************************************
  // Borda pertence a particao de menor numero
  // dentre as particoes de seus elementos
  // ************************************************
  for(int i = 0; i < NBORDER; i++) {
    border[i].part_num = el[border[i].elemento[0]].show_part_num();
  }

  // ****************************************
  // Criar listas de elementos e bordas
  //*****************************************

	// Dimensionar a Mat1<PARTICAO> Particao
	//Particao.resize(NumPart);

  for(int i=0;i<NumPart;++i) {
    PARTICAO Par;
    Par.nele=0;
    Par.nbor=0;
    Par.ngho=0;
    Particao.push_back(Par);
  }

  for(int i=0;i<NELEM;++i) {
    Particao[el[i].show_part_num()].nele++;
  }

  for(int i=0;i<NBORDER;++i) {
    int temp;
    int t=border[i].tipo;
    int k=border[i].part_num;
    Particao[k].nbor++;
    if(t==2) { // borda interior
      temp=el[border[i].elemento[1]].show_part_num();
      if(temp!=k)Particao[k].ngho++;
    }
  }

  for(int i=0;i<NumPart;++i) {
    Particao[i].Aloca_memoria();
    //Particao[i].ele.resize(Particao[i].nele);
    //Particao[i].bor.resize(Particao[i].nbor);
    //Particao[i].gho.resize(Particao[i].ngho);
  }

  int count[NumPart];

  for(int i=0;i<NumPart;++i) count[i]=0;
  for(int i=0;i<NELEM;++i) {
    int k =el[i].show_part_num();
    Particao[k].ele[count[k]++]=i;
  }

  for(int i=0;i<NumPart;++i) count[i]=0;
  int count1[NumPart];
  for(int i=0;i<NumPart;++i) count1[i]=0;
  for(int i=0;i<NBORDER;++i) {
    int temp;
    int t=border[i].tipo;
    int k=border[i].part_num;
    Particao[k].bor[count[k]++]=i;
    if(t==2) { // borda interior
      int ntemp = border[i].elemento[1];
      temp=el[ntemp].show_part_num();
      if(temp!=k)Particao[k].gho[count1[k]++]=ntemp;
    }
  }

  if(myid==0) {
    // Eco da Particao
    FILE * fout_par;
    fout_par=fopen("eco_particao","wb");
    fprintf(fout_par,"NumPart %d\n",NumPart);
    fprintf(fout_par,"Vertices\nindice particao\n");
    for(int i=0;i<NUMNP;++i)fprintf(fout_par,"%d %d\n",i,V[i].part_num);
    fprintf(fout_par,"Elementos\nindice particao\n");
    for(int i=0;i<NELEM;++i)fprintf(fout_par,"%d %d\n",i,el[i].show_part_num());
    fprintf(fout_par,"Bordas\nindice tipo particao\n");
    for(int i=0;i<NBORDER;++i)fprintf(fout_par,"%d %d %d\n",i,border[i].tipo,border[i].part_num);

    // Eco das listas

    fprintf(fout_par,"\nLISTAS\n");
    for(int i=0;i<NumPart;++i) {
      fprintf(fout_par,"Particao %d; nele=%d nbor=%d ngho=%d\n",i,Particao[i].nele,Particao[i].nbor,Particao[i].ngho);
      fprintf(fout_par,"elementos da particao %d\n",i);
      for(int j=0;j<Particao[i].nele;++j) {
        fprintf(fout_par,"ele[%d] = %d\n",j, Particao[i].ele[j]);
      }
      fprintf(fout_par,"elementos ghosts da particao %d\n",i);
      for(int j=0;j<Particao[i].ngho;++j) {
        fprintf(fout_par,"gho[%d] = %d\n",j, Particao[i].gho[j]);
      }
      fprintf(fout_par,"bordas da particao %d\n",i);
      for(int j=0;j<Particao[i].nbor;++j) {
        fprintf(fout_par,"bor[%d] = %d\n",j, Particao[i].bor[j]);
      }
    }
    fclose(fout_par);
  }
};

// *************************************************************************
// *************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::Facet_eco(FILE * fout1)
// *************************************************************************
{
  printf(" template <typename ElemType,int N_VAR,int N_FIELDS> GeProb<ElemType,N_VAR,N_FIELDS>::Facet_eco \n\n");
  fprintf(fout1,"\n******************************************************");
  fprintf(fout1,"**************************");
  fprintf(fout1,"\ntemplate <typename ElemType,int N_VAR,int N_FIELDS> GeProb<ElemType,N_VAR,N_FIELDS>::Facet_eco\n");
  fprintf(fout1,"******************************************************");
  fprintf(fout1,"**************************\n");
  fprintf(fout1,"\nEco das bordas: NBORDER= %4d\n",NBORDER);
  fprintf(fout1,"   i tipo num_elem elem[0] facet[0] elem[1] facet[1]\n");

  for(int i=0; i<NBORDER;++i) {
    fprintf(fout1,"%4d    %d        %d   %4d        %d",i,facet[i].tipo,facet[i].num_elem,facet[i].elemento[0],facet[i].num_local[0]);
    if(facet[i].num_elem==2) fprintf(fout1,"    %4d        %d",facet[i].elemento[1],facet[i].num_local[1]);
    fprintf(fout1,"\n");
  }

  fprintf(fout1,"\n******************************************************");
  fprintf(fout1,"**************************");
  fprintf(fout1,"\ntemplate <typename ElemType,int N_VAR,int N_FIELDS> GeProb<ElemType,N_VAR,N_FIELDS>::Facet_eco\n");
  fprintf(fout1,"******************************************************");
  fprintf(fout1,"**************************\n");
  fprintf(fout1,"\nEco das bordas: NBORDER= %4d\n",NBORDER);
  fprintf(fout1,"   i tipo num_elem elem[0] border[0] elem[1] border[1]\n");
  for(int i=0; i<NBORDER;++i) {
    fprintf(fout1,"%4d    %d        %d   %4d        %d",i,border[i].tipo,border[i].num_elem,border[i].elemento[0],border[i].num_local[0]);
    if(border[i].num_elem==2) fprintf(fout1,"    %4d        %d",border[i].elemento[1],border[i].num_local[1]);
    fprintf(fout1,"\n");
  }
};

// ****************************************************************************
// ****************************************************************************

// *********************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::Ler__Geometria(FILE *finput,double *coord,int *Ed,int *BC)
// *********************************************************************************
{
  // Seq 02.02: DG_Prob::Ler_e_Processar_malha
  // **************************************************************************
  // Ler os dados dos elementos
  // Arquivos no formato .meu diferem do arquivos do Gambit
  // **************************************************************************
  //fscanf(finput,"%d",&NUMNP);
  int n,naux;
  double x,y,z;
  for(int i=0;i<NUMNP;++i) {
    fscanf(finput,"%d %lf %lf %lf",&n, &x, &y, &z);
    if(n>NUMNP || n<0) {
      printf("Erro no arquivo de dados: Vertices\n");
      exit(0);
    }
    else {
      naux=3*n;
      coord[naux++]=x;
      coord[naux++]=y;
      coord[naux  ]=z;
    }
  }
  printf("Leu NUMNP= %d Vertices\n",NUMNP);
  //fscanf(finput,"%d",&NELEM);
  int tipo,numv,v;

  naux=0;
  for(int i=0;i<NELEM;++i) {
    fscanf(finput,"%d %d %d",&n,&tipo,&numv);
    if(n>=NELEM || n<0) {
      printf("Erro no arquivo de dados: Elemento %d\n",n);
      exit(0);
    }
    else {
      //naux=10*n;
      Ed[naux++]=n;
      Ed[naux++]=tipo;
      Ed[naux++]=numv;
      //cout << "echo da entrada de elementos:" << n <<" " <<tipo<< " "<< numv<< endl;
      if(numv>8 || numv<0) {
        printf("Erro no Arquivo de dados: numv %d incompativel\n",numv);
        exit(0);
      }
      else {
        for(int j=0;j<numv;++j) {
          fscanf(finput,"%d",&v);
          Ed[naux++]=v;
        }
      }
    }
  }
  printf("Leu NELEM= %d Elementos\n",NELEM);
  cout<<"Numero de dados lidos em Buffer_E: " << naux-1 <<endl;
  int numel,tipoel,a;
  if(DNBC>0){
    BC[0]=DNBC;
    naux=1;
    for(int i=0;i<DNBC;++i){
      fscanf(finput,"%d %d",&tipo,&n);
      BC[naux++]=tipo;
      BC[naux++]=n;
      for(int j=0;j<n;++j){
        fscanf(finput,"%d %d %d",&numel,&tipoel,&a);
        BC[naux++]=numel;
        BC[naux++]=tipoel;
        BC[naux++]=a;
      }
    }
  }
  else BC[0]=0;
  printf("Leu DNBC= %d Condicoes de contorno\n",DNBC);
  //fclose(finput);
};

// ****************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::MPI_Recebe_Dados_GeProb()
// ****************************************************************************
{
  //MPI::COMM_WORLD.Barrier();
  Comm->Barrier();
  if(comm_size > 1) {

    cout << "Ponto MPI_Recebe_Dados myid " << myid << endl;

    MPI::COMM_WORLD.Bcast(&NUMNP,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&NELEM,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&DNBC,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&NumPart,1,MPI::INT,0);
  }
};
// ****************************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::ResolverComTrilinos(const std::string Package,//const char * Package,
                                                          Epetra_Map Map,
                                                          Teuchos::RCP<Epetra_FECrsMatrix>  A,
                                                          Teuchos::RCP<Epetra_FEVector>  RHS,
                                                          double_t & norm_delta_X)
// ****************************************************************************************
{
  Epetra_MultiVector X(Map,1);
  Epetra_LinearProblem Problem(A.get(), &X, RHS.get());

  if(Package=="Amesos"){
    //printf(" Resolve com Amesos\n");
    // *****************************************
    // Resolucao do sistema Amesos
    // *****************************************
    Amesos_BaseSolver * Solver;
    Amesos Factory;
    // std::string AmesosSolverType = "Amesos_Umfpack";
    //"Amesos_Superludist";
    //"Amesos_Umfpack" ;

    bool IsAvailable = Factory.Query(AmesosSolverType);
    assert(IsAvailable);

    Solver = Factory.Create(AmesosSolverType,Problem);
    //MPI::COMM_WORLD.Barrier();
    Comm->Barrier();
    Solver->SymbolicFactorization();
    Solver->NumericFactorization();
    Solver->Solve();
    delete Solver; Solver=nullptr;

    // **************************************
    // Fim da Resolucao Amesos
    // **************************************
  }

  else if (Package=="AztecOO"){
    //printf(" Resolve com AztecOO\n");
    // *****************************************
    // Resolucao do sistema AztecOO
    // *****************************************
    AztecOO solver (Problem);
    solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    solver.SetAztecOption(AZ_output, AZ_none);
    solver.SetAztecOption(AZ_solver, AZ_gmres);
    solver.Iterate(1000,1.0e-08);
    std::cout << "Solver performed " << solver.NumIters() << " iterations." << std::endl
              << "Norm of true residual = " << solver.TrueResidual() << std::endl;
    // **************************************
    // Fim da Resolucao AztecOO
    // **************************************
  }

  else {
#ifdef HAVE_MPI
    MPI::Finalize();
#endif
    exit(0);
  }
  X.NormInf(&norm_delta_X);
  // ***********************
  // Exemplo da difusao de X0 (=copia de X)
  // processador 0 recebe as partes dos demais processadores
  // monta X0 e faz Broadcast de seu conteudo completo
  std::vector<double>  X0  (NumD); //hoje
#ifdef HAVE_MPI
  int NumMyVars_target;
  if(myid==0)
    NumMyVars_target=NumD;
  else
    NumMyVars_target=0;
  Epetra_Map TargetMap(-1,NumMyVars_target,0,*Comm);
    Epetra_Export Exporter(Map,TargetMap); //analogo a Montagem da matriz global; TargetMap eh 1-1, Map  pode nao ser
     //Epetra_Import Importer(TargetMap,Map); //analogo a distribuicao da solucao para os elementos Map eh 1-1
                                            // TargetMap pode nao ser
  Epetra_Vector Y(TargetMap);

    Y.PutScalar(0.0);
  Y.Export(X,Exporter,Add);
  //Y.Import(X,Importer,Insert);
    Y.ExtractCopy(&X0[0]);

  MPI::COMM_WORLD.Bcast(&X0[0],NumD,MPI::DOUBLE,0);
#else
  for ( int i = 0; i < NumD; i++) {
    X0[i]=X[i];
  }
#endif

  // ***************************************
  // Avanca a solucao
  // ***************************************
  for ( int i = 0; i < NELEM; i++)
    el[i].Avancar_u0(&X0[0]);
  //delete [] X0; X0=nullptr;
};

// ***************************************************************************
// ***************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::Ler_e_Processar_malha(char *arq_geo)
{
  //Seq 01.02: GeProb::Ler_e_Processar_malha
  // *************************************************************************
  // Ler os parametros para criar a  malha do problema:  GeProb
  // *************************************************************************
  FILE * finput_geo, *finput_part;
  char arq_part[256]; // arquivo da particao da malha

  if(myid==0) {
    strcpy(arq_part,arq_geo); // copia o nome arq_geo em arq_part
    strcat(arq_part,"_part\0"); // nome do arquivo contendo a particao tem terminacao _part

    finput_geo=fopen(arq_geo,"rb"); // Arquivo de geometria

    fscanf(finput_geo,"%d %d %d",&NUMNP,&NELEM,&DNBC);


    if((finput_part=fopen(arq_part,"rb"))!=NULL) {
      printf("Lendo arquivo de Particao\n");

      fscanf(finput_part,"%d",&NumPart);
    }
    else {
      printf("Nao existe arquivo de Particao %s!!!!\n",arq_part);
      NumPart=1;
    }

    printf("myid %d Parametros NUMNP = %d NELEM = %d DNBC = %d NumPart = %d\n",myid,NUMNP,NELEM,DNBC,NumPart);
  } // myid == 0

#ifdef HAVE_MPI
  // *************************************************************************
  // Broadcast os parametros da malha geometrica e dos espacos de funcoes
  // do problema
  // *************************************************************************
  // MPI::COMM_WORLD.Barrier();
  Comm->Barrier();
  if(comm_size > 1) {
    cout << "Ponto MPI_Recebe_Dados myid " << myid << endl;

    MPI::COMM_WORLD.Bcast(&NUMNP,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&NELEM,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&DNBC,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&NumPart,1,MPI::INT,0);
  }
  Comm->Barrier();
#endif

  V = new Vertice [NUMNP];
  el= new ElemType [NELEM];

  // ********************************************************
  // Se o NumPart for incompativel com o numero de processos
  // Parar a execucao
  // ********************************************************
  if(comm_size != NumPart) { // Nao faz a iteracao
    if(myid==0) {
      printf("\nGeProb<>::Ler_e_Processar_malha(): Numero de processadores eh diferente do numero de particoes\n");
      printf("\nTerminar sem nada fazer\n\n");
    }

#ifdef HAVE_MPI
    MPI::Finalize();
# endif

    exit(0);
  } // comm_size != NumPart
  // ********************************************************

  // *************************************************************************
  // Ler os dados da  malha do problema
  // *************************************************************************
  int controle[4];
  controle[0]=3*NUMNP;
  controle[1]=11*NELEM;
  controle[2]=122*DNBC;
  controle[3]=NumPart + NUMNP;

  /*
    cout << "controle[0] = "<< controle[0] << endl;
    cout << "controle[1] = "<< controle[1] << endl;
    cout << "controle[2] = "<< controle[2] << endl;
    cout << "controle[3] = "<< controle[3] << endl;
  */

  double buffer_V [controle[0]];// maior que 3*NUMNP
  int    buffer_E [controle[1]];// maior que 11*NELEM
  int    buffer_BC[controle[2]];// maior que 10*DNBC
  int    buffer_Pa[controle[3]];// Numero de dados da particao +1
  // ********************************************************

  if(myid==0) {

    Ler__Geometria(finput_geo,buffer_V,buffer_E,buffer_BC); //
    fclose(finput_geo);

    if(finput_part!=NULL) { // **********  Ler arquivo com particao *******
      for(int i=0;i<controle[3];++i) {
        fscanf(finput_part,"%d",&buffer_Pa[i]);
      }
      fclose(finput_part);
    }
  } // if(myid==0)
  // ********************************************************

#ifdef HAVE_MPI
  // *************************************************************************
  // Broadcast dados geometricos da malha do problema
  // *************************************************************************

  if(comm_size > 1) {
    // MPI::COMM_WORLD.Barrier();
    Comm->Barrier();
    MPI::COMM_WORLD.Bcast(buffer_V, controle[0],MPI::DOUBLE,0); // Vertices
    // MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Bcast(buffer_E, controle[1],MPI::INT,0); // elementos
    // MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Bcast(buffer_BC,controle[2],MPI::INT,0); // Condicoes de contorno
    // MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Bcast(buffer_Pa,controle[3],MPI::INT,0); // Particao
  }
#endif

  // *************************************************************************
  // Terminou leitura e propagacao dos dados
  // *************************************************************************

  // *************************************************************************
  // Começar a processar os dados de entrada
  // *************************************************************************

  // *************************************************************************
  // Criar os vertices e armazenar as coordenadas
  // *************************************************************************

  int n=0;
  for(int i=0;i<NUMNP;++i) { // *** Coordenadas dos Vertices  ****************
    V[i].x=buffer_V[n++];
    V[i].y=buffer_V[n++];
    V[i].z=buffer_V[n++];
  }
  if(myid==0) printf("Releu o NUMNP= %d\n",NUMNP);

  // *****************************************************
  // ***** Construir Elementos Fisicos (PhElem) **********
  // *****************************************************
  // 10/07/2014
  std::vector< std::vector<int> > face_mask_vec;
  std::vector<int> face_mask_local (5,0);

  int ii=0;
  for(int i=0;i<NELEM;++i) {
    int ind,tipo,numv,vert[8];
    ind =buffer_E[ii++];//i;
    tipo=buffer_E[ii++];
    numv=buffer_E[ii++];
    // cout << "numv teste:"<<numv<<endl;
    for(int j=0;j<numv;++j) {
      vert[j]=buffer_E[ii++];
    }
    // el[i].set_NumLocalVars(1);// (N_FIELDS);
    el[i].set_ptvert(&V[0]);
    el[i].set_type(tipo);
    if(tipo==4) {
      tetrahedro_faz_face_mask(numv,vert,face_mask_local);
      face_mask_local[4]=i; // numero do elemento a quem pertence
      face_mask_vec.push_back(face_mask_local);
    }
    if(tipo==5) { // hexaedro requer ordenacao dos nos lidos de gmsh;
        int temp = vert[2];
        vert[2]=vert[3];
        vert[3]=temp;
        temp   = vert[6];
        vert[6]= vert[7];
        vert[7]= temp;
      }
    el[i].set_Vert_map(numv,vert);
  }

  if(myid==0) printf("Releu o NELEM= %d\n",NELEM);

  // ----------------------FIM-----------------------------
  // ***** Construir Elementos Fisicos (PhElem) **********
  // *****************************************************

  // *****************************************************
  // Construcao das aproximacoes das funcoes
  // *****************************************************
  // processar dados da malha
  // *****************************************************
  // cout << "N_var lido: "<< N_var<< ' '<< FieldOfVar[0]<< ' ' << FieldOfVar[1] << '\n';


  Processar_elementos();
    cout << "Passou Processar_elementos();\n";
  Construir_bordas();
    cout << "Passou Construir_bordas();\n";
  Particionar_malha( buffer_Pa );
  Marcar_condicoes_contorno( buffer_BC,face_mask_vec );

  cout << "dimensao em Ler_e_processar_malha "<< dim << endl;
};

// *****************************************************************************
// 08/07/2014 \/ \/ \/ \/ \/ \/ \/ \/
// *****************************************************************************
// Calcula os valores iniciais para numerar os modos sobre os elementos da malha
// *****************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::mesh_suporte(const int & P, int Temp_V[], int Temp_A[], int Temp_F[],int & count)
{
  // Vertices
  for(int i=0;i<NUMNP; ++i)
    Temp_V[i] = count++;

  // Arestas
  for(int i=0;i<NL; ++i) {
    int inc = (P-1);
    if(inc > 0) {
      Temp_A[i] = count;
      count += inc;
    }
    else Temp_A[i]=-999; // sem numeracao
  }
  // Faces
  for(int i=0; i<NF; ++i){
    int inc=0;
    if(Face[i]._tipo == 2) inc = (P-1)*(P-2)/2; // triangulo
    if(Face[i]._tipo == 3) inc = (P-1)*(P-1);// quadrilatero
    if(inc > 0) {
      Temp_F[i] = count;
      count += inc;
    }
    else Temp_F[i]=-999; // sem numeracao
  }
  cout << "mesh_suporte (NUMNP,NL,NF) = (" << NUMNP << ", "<< NL <<", "<< NF <<")" << endl;
  cout << "count = "<< count<< endl;
};

// *****************************************************************************
// Mapeamento global para funcoes continuas atraves das bordas
// Caso Geral. Funciona para todas as dimensoes
// *****************************************************************************

template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::gbnmap_continuous(int & count)
{
  for(int i=0;i<N_VAR;++i){
    gbnmap_continuous(i,count); // overload para cada variavel
  }
  // verificar se interfere no resto do programa
  NumD=count;
};

// *****************************************************************************
// Mapeamento global para funcoes continuas atraves das bordas
// Especializado para a variavel ivar
// Caso Geral. Funciona para todas as dimensoes
// *****************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::gbnmap_continuous(const int & ivar,int & count)
{
  // aloca memoria de vec_map ( mapeamento dos vertices)
  int vert_vec[NUMNP], aresta_vec[NL],face_vec[NF];
  int P = el[0].show_ptr_stdel(ivar)->P_val(0);
  //cout << "contador de modos inicial = "<< count << endl;
  mesh_suporte(P, vert_vec, aresta_vec, face_vec, count);
  //cout << "contador de modos depois de mesh_suporte = "<< count << endl;
  //cout << "dimensao = "<< dim << endl;
  printf("vert_vec de var %d\n",ivar);
  for (int k=0;k<NUMNP; ++k)
    printf("%d %d\n",k,vert_vec[k]);

  printf("\naresta_vec de var %d\n",ivar);
  for (int k=0;k<NL; ++k)
    printf("Aresta %d modo inicial = %d Na=%d Nb=%d\n",k,aresta_vec[k], Aresta[k].Na, Aresta[k].Nb);

  printf("\nface_vec de var %d\n",ivar);
  for (int k=0;k<NF; ++k)
    printf("%d %d\n",k,face_vec[k]);

  for(int e=0;e<NELEM;++e) el[e].gbnmap_vertices(vert_vec,ivar);

  if(dim > 1) {
    for(int e=0;e<NELEM;++e) el[e].gbnmap_arestas(aresta_vec,ivar);
    if(dim == 3)
      for(int e=0;e<NELEM;++e) el[e].gbnmap_faces(face_vec,ivar);
  }

  // construir mapeamento interior
  for(int e=0;e<NELEM;++e){
    // cout << "elemento "<< e<< " count = "<< count << ".... ";
    el[e].gbnmap_interior(ivar,count);
    //cout << count << endl;

  }
  //cout << "contador de modos no final da variavel "<< ivar <<  " = " << count << endl;
};
// ****************************************************************************
// Mapeamento global para funcoes descontinuas atraves das bordas
// Especializado para a variavel ivar
// ****************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::gbnmap_discontinuous(const int & ivar,int & count)
{
  for(int i=0;i<NELEM;++i)
    el[i].inicia_gbnmap(ivar,count);
};

// *******************************************************************
// Funcao requer muita atencao na sua elaboracao
// Constroi as bordas a partir dos dados contidos nos elementos
// fisicos PhElem el[]
// *******************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::Construir_bordas()
{
  // Referencia cruzada elemento <-> border
  cout << "template <typename ElemType,int N_VAR,int N_FIELDS> GeProb<ElemType,N_VAR,N_FIELDS>::Construir_bordas\n";
  // Dimensao espacial do problema é dado pela dimensao do primeiro elemento fisico.
  dim = el[0].show_ptr_stdel(0)->ndim_val();
  cout << "dimensao em template <typename ElemType,int N_VAR,int N_FIELDS> GeProb<ElemType,N_VAR,N_FIELDS>::Construir_bordas "<<dim << endl;
  switch(dim){

    case 1:

      NBORDER = NUMNP;
      // Iniciar os borders
      border = new EDGE [NBORDER];

      for(int i=0;i<NBORDER;++i){
        //EDGE borda;
        border[i].tipo=0;
        border[i].num_elem=0;
        border[i].elemento[0]=0;
        border[i].elemento[1]=0;
        border[i].pdir=nullptr;
        border[i].sdir=nullptr;
        //border.push_back(borda);
      }

      for(int i=0; i<NELEM; ++i) {
        int numv = el[i].show_numv();
        for(int j=0;j<numv;++j){
          int n_border = el[i].show_Vert_map(j);
          int n = border[n_border].num_elem;
          border[n_border].elemento[n]=i;
          border[n_border].num_local[n]=j;
          border[n_border].num_elem++;
          el[i].set_border_num(j,n_border);//19/07/2014
        }
      }

      break; // fim de dimensao 1

    case 2:

      // Testado com quadrilateros
      cout << "Caso 2 (dim) "<< endl;
      NBORDER = NL;
      border = new EDGE [NBORDER];
      // Iniciar os borders
      for(int i=0;i<NBORDER;++i){
        //EDGE borda;
        border[i].tipo=0;
        border[i].num_elem=0;
        border[i].elemento[0]=0;
        border[i].elemento[1]=0;
        border[i].pdir=nullptr;
        border[i].sdir=nullptr;
        //border.push_back(borda);
      }

      for(int i=0; i<NELEM; ++i) {
        int nume = el[i].show_nume();
        for(int j=0;j<nume;++j){

          int n_border = el[i].show_Aresta(j);
          int n = border[n_border].num_elem;

          border[n_border].elemento[n] =i;
          border[n_border].num_local[n]=j;
          border[n_border].sinal[n]    =el[i].show_sinal(j);
          if(n==0) { // primeira vez que a borda e visitada
            int n0 = el[i].show_ptr_stdel(0)->aresta_lvert(j,0);//aresta j vertice 0
            int n1 = el[i].show_ptr_stdel(0)->aresta_lvert(j,1);//aresta j vertice 1
            int V0 = el[i].show_Vert_map(n0); // numero global do vertice 0
            int V1 = el[i].show_Vert_map(n1); // numero global do vertice 1
            if (V0 > V1) { // borda e definida sempre do menor vertice para o maior
              int aux = V0;
              V0 = V1;
              V1 = aux;
            }
            border[n_border].Na=V0;
            border[n_border].Nb=V1;
            // calcular comprimento, w e Jacobiano
          }
          border[n_border].num_elem++;
          el[i].set_border_num(j,n_border);//19/07/2014
        } // for(j =0; j<nume)
      }

      break; // fim de dimensao 2

    case 3:

      NBORDER = NF;
      // Iniciar os borders
      border = new EDGE [NBORDER];

      for(int i=0;i<NBORDER;++i){
        //EDGE borda;
        border[i].tipo=0;
        border[i].num_elem=0;
        border[i].elemento[0]=0;
        border[i].elemento[1]=0;
        border[i].pdir=nullptr;
        border[i].sdir=nullptr;
        //border.push_back(borda);
      }

    
      for(int i=0; i<NELEM; ++i) {
        int numf = el[i].show_numf();
        for(int j=0;j<numf;++j){
          int n_border = el[i].show_Face(j);
          int n = border[n_border].num_elem;
          border[n_border].elemento[n]=i;
          border[n_border].num_local[n]=j;
          border[n_border].num_elem++;
          el[i].set_border_num(j,n_border);//19/07/2014
          // calcular area, w e Jacobiano
        }

      }

      break;

    default:
      printf("Elemento lido nao eh de tipo conhecido\n");
      exit(0);
  }

  // determinar o tipo (default = 0)
  for(int i=0;i<NBORDER;++i){
    if(border[i].num_elem == 2)
      border[i].tipo=2;
  }

  // ***************************************************
  // Construir vetor normal aa superficie das bordas
  // ***************************************************
  for(int i=0;i<NBORDER;++i){
    double normal[3];
    double area;
    int num_local = border[i].num_local[0];
    int e = border[i].elemento[0];
    el[e].vetor_superficie(num_local,area,normal);
    border[i].comprimento=area;
    for(int k=0;k<dim;++k) {
      border[i].normal[k]=normal[k];
    }
  }

};

// *****************************************************
// Marcar as condicoes de contorno no GeProb
// *****************************************************
// ****************************************************************************************
template <typename ElemType,int N_VAR,int N_FIELDS>
void GeProb<ElemType,N_VAR,N_FIELDS>::Marcar_condicoes_contorno(int *BC,
                                                                std::vector< std::vector<int> > face_mask)
// ****************************************************************************************
{
  // Inicializar o bflag com valores 0
  for(int i=0;i<NG;++i) bflag.push_back(0); //bflag[i]=0;//bflag=0: conhecido

  int naux,tipo,n;
  int elnum,eltype,facenum;//newfacenum;

  DNBC=BC[0];
  naux=1;

  for(int k=0;k<DNBC;++k){// BOUNDARY CONDITIONS
    tipo=BC[naux++]; // tipo de condicao de contorno
    n=BC[naux++]; //< numero de elementos que tem essa condicao de contorno

    for(int j=0;j<n;++j){
      elnum=BC[naux++]; //< numero do elemento
      eltype=BC[naux++]; //< tipo do elemento
      facenum=BC[naux++]; //< face do elemento onde vale a condicao de contorno

      // caso especial para o tetraedro
      // Nao testado ainda. 01/11/2016
      if(eltype==4){
        int flag=0;
        for(int i=0;i< face_mask.size() && flag==0; ++i) {
          if(elnum==face_mask[i][4]) {
            facenum=face_mask[i][facenum]; // Tetraedro
            flag=1;
          }
        }
      }

      el[elnum].set_border_tipo(border,facenum,tipo); //incluido 01/nov/2016
    } // incluido 01/nov/2016
  } // for(int k=0; k< DNBC; ++k) incluido 01/nov/2016
};
// ******************************************************************************************

#endif /* _GeProb_headers */
// ****************************************************


