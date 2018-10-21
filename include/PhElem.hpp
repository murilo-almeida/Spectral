//
//  PhElem.hpp
//  SDG
//  Classe base
//  Elementos fisicos PhElem definido por template
//  para poder especificar o numero de variaveis no momento
//  da sua definicao.
//
//  Created by Murilo Almeida on 24/10/16.
//  Copyright © 2016 Murilo Almeida. All rights reserved.
//

#ifndef PhElem_hpp
#define PhElem_hpp

#include "spectral.h"

// *************************************************************************
template <int NumVariaveis  >
class PhElem
{
public:

  PhElem(); //(const int n = 1);
  ~PhElem(){/*cout << "destruir PhElem\n";*/ };
  void inicia_vetores();
  void finaliza_vetores();
  int show_NumLocalVars(){return(NumLocalVars);};;
  void teste();
  void set_ptr_stdel(Stdel * pointer, Stdel * pointer1);
  void set_ptr_stdel(Stdel * pointers[]);
  void set_ptr_stdel(Stdel * pointer);
  void set_ptr_stdel_var(const int ind, Stdel * pointer);
  void set_Vert_map(const int & n_in,int ver_temp[]);
  void set_ptvert(const Vertice * pointer);
const Vertice * show_ptvert(){return ptvert;};
  void read_vertices(FILE *, const int & );
  void set_sgn();
  int show_sgn(const int & ia, const int & i){return (sgn[ia][i]);};
  const int show_gbnmap(int ivar, int modo) const { return gbnmap[ivar][modo];};
  int show_numv(){return(numv);};
  int show_nume(){return(nume);};
  int show_numf(){return(numf);};
  void set_type(const int & num);
  int type_val();
  void print_numeracao(FILE * fout,const int & k);
  void print_modes(FILE * fout,const int & k) const;
 // void inicia_vetores();
  //void inicia_tracos(EDGE *  border);
  //void inicia_funcoes_na_borda( EDGE *  border);
  //void finaliza_vetores();
  void zera_vetores(const int & k);
  void zera_bs(const int & k);
  void make_vector(const int & ia,double (*f)(double,double,double));
  void make_vector_Elast2D_dyn(const double a[]);
  void make_Kcomp_Elast2D(const double lambda, const double mu);
  void make_Kcomp_Elast2D_dyn(const double lambda,const double mu,
                              const double a[]);
  void make_K_Elast2D(const double lambda, const double mu,double ** K);
  void VectorElast2D(const int & i,double sum[]);
  void VectorElast2D(double vector[]);
  //void AlocarMatrizesK(Matrix * ptKcomp,Matrix * ptKi_inv,Matrix *ptKcKi_inv);
  //void AlocarMatrizesM(Matrix * pt1,Matrix * pt2,Matrix *pt3);
  void P_eval_print_Elast2D(FILE * fout,
                            const double X[],double lambda, double mu);
  void P_eval_print_Elast2D_dyn(FILE * fout,const int & prnflag,
                                const double X[],double lambda, double mu,
                                const double a[]);
  void make_K_Imiscivel2F(const double,const double,const double, const double,
                          const double [],double (*)(double),
                          double (*)(double),double (*)(double));
  double Kcomp(const int & i, const int & j);
  double K_Imisc2F(const int & i, const int & j);
  void make_Vector_Eq_S1(const double, const double,
                         const double, const double,
                         const double [3],double (*)(double));
  //int map(const int & i);
  int map(const int & var,const int & no_local);
  void mapa_inverso(const int & ivar, const double X[]);
  //double MassMatrix_value(int i,int j);
  void P_eval_u(const double X[]);
  void P_eval_u(const double X[],const int & ia);
  void P_eval_u_p1(const double X[],const int & ia);
  void P_eval_print(const double X[],const int & numf,FILE *fout,
                    double (*)(double,double,double));
  void P_print(const int & numf, FILE *file);
  void P_eval_vert(const double X[],double f_vert[]);
  void P_eval_phys(const int & numf,double f0[]);
  void Increment_field(const double X[],const int & numf, const double dt);
  void compute_JV(const int & ia);
  void print_matrices(FILE * fout);
  /*
   void Processar_dados(const int & nel,int& NG,int& NL,std::vector<EDGE>& border,int Ng[],
   int & NF, int Face[], int Fng[], int f_mask[]);
   */
  void Processar_dados(int& NL, std::vector<ARESTA>& aresta,
                       int & NF, std::vector<FACE> & face_vec);
  void check_connectivity(FILE * fout,const int & ia);
  void assign_gbnmap(const int & ia,const int & i, const int & val);
  void inicia_gbnmap(int & count);
  void inicia_gbnmap(const int & ivar,int & count);
  void inicia_gbtrbmap(int & count);
  void set_gbnmap(const int & ia,const int gbnmap_temp[],
                  const int sgn_temp[]);
  void set_stgbtrbmap(const int & b, const int & vmapM, const int & vmapP);
  int show_gbtrbmap(){ return(gbtrbmap); };
  int show_stgbtrbmapM(const int & a){ return(stgbtrbmapM[a]); };
  int show_stgbtrbmapP(const int & a){ return(stgbtrbmapP[a]); };
  void check_gradiente(FILE * fout,const int & ia);
  Stdel * show_ptr_stdel(const int & n) const;

  //void make_bflag(int bflag[], int facenum);
  void BoundCond(const int & face, int bflag[], double bc[],
                 double (*func)(double,double,double),
                 const int & nf, const int & bndtype);
  void Neumann(const int & aresta,const int & varn,
               const double rho1,const double mu1,
               const double rho2,const double mu2,
               const double gravidade[]);
  double show_u(const int & i);
  double show_u(const int &,const int &);
  double show_b0(const int & i);
  double show_b0(const int &,const int &);
  double show_bs(const int &,const int &);
  //void set_mass_density(double);
  //void set_porosidade(double);
  //void set_fontes(double sw,double sn){qw=sw; qn=sn;};

  //double show_rho(){return rho;};
  double show_J(){return J;};
  //void set_permeabilidade(const double,const double, const double);
  void make_MassMatrices();
  int show_Vert_map(int i){return Vert_map[i];};
  int show_Aresta(int i){return aresta_map[i];};
  int show_Face(int i){return face_map[i];};
  int show_sinal(int i){return sinal[i];};
  int show_border_num(int i){return border_num[i];};
  int show_part_num(){return part_num;};
  //double show_Volume(){return Volume;};
  void set_part_num(const int & num = -1);
  double Phi_val(const int & var, const int & ind, const int & pos){
    return ptr_stdel[var]->show_Phi_val(ind,pos);};
  void projetar_C0(FILE *file,double (*func)(double,double,double),
                   const int & ivar);
  void transformacao_direta(double f[],const int & ivar);
  void teste_transformacao_direta(FILE * fin, FILE * fout, const int & npoints,const double coord[]);

  void set_border_tipo( EDGE * border,const int & aresta,const int & t);
  void set_border_num(const int  & aresta, const int & num);
  //void Atualizar_valores(FILE * fout = NULL);// NULL eh o valor default
  void Imprimir_valores(FILE * fileout, const int & npoints, const double coord[]);
  void printwGQJ(FILE * fileout);
  void Avancar_u0(const double X[],const double relax = 1.0);
  void Atualizar_u0(const double X[]);
  void Copia_u0_em_(double X[]);
  void Copia_u0_em_(Teuchos::RCP<Epetra_Vector> X);
  void Comparar_u0(const double X[]);
  void Salvar_u0();
  void Restaurar_u0();
  void escrever_restart(FILE * fout);
  void ler_restart(FILE * fin);
  void ler_restart_buffer(FILE * fin,double * buff, int & count);
  void restart_element(double * buff,int & count);
  //void testar_traco(FILE * f_eco);
  //void echo_traco(FILE * f_eco = nullptr);
  void teste_gradiente();
  // void teste_tracoI(Fluids fls);
  int  get_trace_border_map(const int & aresta);
  // ***********************************************
  void gbnmap_vertices(std::vector< std::vector<int>> gbn);//int ** gbn);
  void gbnmap_vertices(int * gbn,const int & ivar);
  void gbnmap_faces(const int face_vec[],const int & ivar);
  void gbnmap_aresta(const int & aresta,const int & sinal,const int & ivar,
                     const int & count,int &inc);
  void gbnmap_arestas(const int aresta_vec[],const int & ivar);// nova implementacao
  void gbnmap_interior(const int & ivar,int &count);
  void anexa_gbnmap(const int & ivar, vector <int> &);

  void vetor_superficie(const int & num_local,double & area, double normal[3]);


  //private:
protected:

  Stdel * ptr_stdel[NumVariaveis];
  int type;
  const Vertice * ptvert;
  int NumLocalVars = NumVariaveis;
  int ndim; //!< Spatial dimension
  int numv; //!< Number of vertices
  int nume; //!< Number of edges
  int numf; //!< Number of faces
  int numborders; //!< Number of borders in the element
  int numn[NumVariaveis];
  int numb[NumVariaveis];
  int Vert_map[8], aresta_map[12],face_map[6], border_num[12], sinal[12];//!< sinal refers to the borders
  int part_num; //!< numero da particao a qual pertence
  int gbnmap[NumVariaveis][MAXNN]; //!< mapping from local to global nodes (or modes)
  int sgn[NumVariaveis][MAXNN];  //!< sign of the local mode (or mode)
  double J; //!< Jacobian from the physical to the standard element
  double * JV; //!<Ponteiro para matriz do Jacobiano nos pontos gaussianos;
  double * b0[NumVariaveis];     //!<[MAXMODES*MAXNFIELDS];  vector of body forces
  double * bs[NumVariaveis];     //!< vector of force boundary conditions
  double * u0[NumVariaveis]; //!<vector with the coefficients of each mode
  double * usave[NumVariaveis];
  //double * ua[NumVariaveis];

  int gbtrbmap ;/*!< \brief Beginning index on the global array where
                 * the local border Gauss quadrature points (traces)
                 * are mapped; this map facilitates
                 * the use of mpi for paralell calculations over
                 * element borders */
  int * stgbtrbmapM; /*!< \brief Array containing the start point in the global
                      *  trace array of the internal trace of the border */
  int * stgbtrbmapP; /*!< \brief Array containing the start point in the global
                      *  trace array of the external trace of the border */
  int vetores_iniciados; // = 1; indica que os vetores locais foram iniciados e necessitam ser finalizados

};
/*! \class PhElem
 * Physical Elements
 */


// ****************************************************************************
// Class PhElem ---  Para escoamentos de 2 fluidos imisciveis DG
// ****************************************************************************

template < int NumVariaveis >
PhElem<NumVariaveis>::PhElem()
{
  //set_ptr_stdel(NULL,NULL);
    vetores_iniciados = 0;
  // qn=0.0;
  // qw=0.0;
  // default value
  //printf("Inicializa PhElem: NumLocalVars = %d rho = %g\n",n,rho);
};
// ****************************************************************************

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::escrever_restart(FILE * fout)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=numn[ivar];
    for(int i=0;i<nn;i++) fprintf(fout," %e",u0[ivar][i]);
  }
  fprintf(fout,"\n");
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::ler_restart(FILE * filein)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=numn[ivar];
    for(int i=0;i<nn;i++) fscanf(filein,"%lf",&u0[ivar][i]) ;
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::ler_restart_buffer(FILE * filein,double * buff, int & contador)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=numn[ivar];
    for(int i=0;i<nn;i++) fscanf(filein,"%lf",&buff[contador++]) ;
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::restart_element(double * buff,int & conta)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=numn[ivar];
    for(int i=0;i<nn;i++) u0[ivar][i] = buff[conta++] ;
  }
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::mapa_inverso(const int & ivar, const double X[])
{
  double temp;
  int nn=numn[ivar];
  for(int i=0;i<nn;i++){
    temp = X[gbnmap[ivar][i]];
    u0[ivar][i]= temp;
    if(std::isnan(temp))cout << "Mapa inverso Float was Not a Number: ivar " << ivar << " modo "<< i << endl;
  }

};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Copia_u0_em_(double X[])
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    for(int i=0;i<numn[ivar];i++) X[gbnmap[ivar][i]] = u0[ivar][i];
  }
};
template<int NumVariaveis>
void PhElem<NumVariaveis>::Copia_u0_em_(Teuchos::RCP<Epetra_Vector> X)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    for(int i=0;i<numn[ivar];i++) (*X)[gbnmap[ivar][i]] = u0[ivar][i];
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::anexa_gbnmap(const int & ivar, vector <int> & list)
{
  for(int i=0;i<numn[ivar];i++) list.push_back( gbnmap[ivar][i] );
};


// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::printwGQJ(FILE * fileout)
{
  if(fileout != NULL){
    ptr_stdel[0]->printwGQtofile(fileout,ptvert,Vert_map,JV);
    fflush(fileout);
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Salvar_u0()
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    for(int i=0;i<numn[ivar];++i) usave[ivar][i] = u0[ivar][i];
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Restaurar_u0()
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    for(int i=0;i<numn[ivar];i++) u0[ivar][i] = usave[ivar][i];
  }
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Avancar_u0(const double X[],const double relax)
{
  // default relax = 1.0

  // atualizar os coeficientes: u0 -= X

  for(int ivar=0;ivar<NumVariaveis;ivar++){
    // *****************************************************************
    for(int i=0;i<numn[ivar];i++){                                          // *
      u0[ivar][i] -= relax * X[gbnmap[ivar][i]];                    // *
    }                                                               // *
    // *****************************************************************
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Comparar_u0(const double X[])
{
  // default relax = 1.0



  // comparar os coeficientes: u0 e X

  for(int ivar=0;ivar<NumVariaveis;ivar++){
    // *****************************************************************
    for(int i=0;i<numn[ivar];i++){                                          // *
      if(u0[ivar][i] != X[gbnmap[ivar][i]])
        cout << "Falhou em VarGlobal " << gbnmap[ivar][i] << "\n" ; // *
    }                                                               // *
    // *****************************************************************
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Atualizar_u0(const double X[])
{

  // atualizar os coeficientes: u0 = X

  for(int ivar=0;ivar<NumVariaveis;ivar++){
    // *****************************************************************
    for(int i=0;i<numn[ivar];i++){                                          // *
      u0[ivar][i] = X[gbnmap[ivar][i]];                             // *
    }                                                               // *
    // *****************************************************************
  }
};

/*
 // ****************************************************************************
 template<int NumVariaveis>
 void  PhElem::testar_traco(FILE * f_eco)
 {
 const int qmax=ptr_stdel[0]->qborder_val();
 const int ndim=ptr_stdel[0]->ndim_val();
 int NGQP=ptr_stdel[0]->NGQP_val();
 double * gphi_r[ndim];
 double * gphi_l[ndim];
	double A[2][qmax];
	double B[2][qmax];
 for(int i=0; i<2; i++){
 gphi_r[i]=A[i];
 gphi_l[i]=B[i];
 }

 fprintf(f_eco,"Phi para o elemento\n");
 for(int i=0;i<2;i++){
 fprintf(f_eco,"\n variavel %d\n",i);
 int nn=ptr_stdel[i]->nn_val();
 int nborder=ptr_stdel[i]->nborder_val();
 for(int j=0;j<nn;j++){
 fprintf(f_eco,"         modo %d\n",j);
 fprintf(f_eco,"Grad_Phi[%d][%d]\n",i,j);
 for(int m=0;m<NGQP;m++){
 fprintf(f_eco," m= %d ",m);
 for (int ndir=0; ndir< ndim; ndir++){
 fprintf(f_eco,"%11.4e ",GradPhi[i][j][ndir][m]);
 }
 fprintf(f_eco,"\n");
 }
 printf("AQUI1 *********************************************************\n");
 for(int h=0;h<nborder;h++){
 Traco_grad_phi(h,i,j,gphi_r);
 fprintf(f_eco,"  Traco na aresta %d \n",h);
 for(int m=0;m<qmax;m++){
 for (int ndir=0; ndir< ndim; ndir++){
 fprintf(f_eco,"%11.4e ",gphi_r[ndir][m]);
 }
 fprintf(f_eco,"\n");
 }
 }
 printf("AQUI2 *******************************************************\n");
 }
 }
 // **********************************************************
 // Fim do echo dos resultados - testar_traco
 // **********************************************************
 };
 */


// ****************************************************************************
//template<int NumVariaveis>
//void PhElem<NumVariaveis>::set_NumLocalVars(const int & n)
//{
// NumLocalVars=n;
//};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel(Stdel * point0 )
{
  ptr_stdel[0]=point0;
  numn[0]=ptr_stdel[0]->nn_val();
  numb[0]=ptr_stdel[0]->nb_val();
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel(Stdel * point0, Stdel * point1 )
{
  if (NumVariaveis > 2) {
    ptr_stdel[0]=point0;
    ptr_stdel[1]=point1;
    numn[0]=ptr_stdel[0]->nn_val();
    numb[0]=ptr_stdel[0]->nb_val();
    numn[1]=ptr_stdel[1]->nn_val();
    numb[1]=ptr_stdel[1]->nb_val();

  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel_var(const int var, Stdel * point0 )
{
  ptr_stdel[var]=point0;
  numn[var]=ptr_stdel[var]->nn_val();
  numb[var]=ptr_stdel[var]->nb_val();
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel(Stdel *  ponteiros[])
{
  for(int i=0;i<NumVariaveis;++i){
    ptr_stdel[i]=ponteiros[i];
    numn[i]=ptr_stdel[i]->nn_val();
    numb[i]=ptr_stdel[i]->nb_val();
  }
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptvert(const Vertice * pointer)
{
  ptvert=pointer;
};


// ****************************************************************************
// Especifica o tipo(label) da borda
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_border_tipo(EDGE * border,
                                           const int & aresta,const int & t)
{
  int en = border_num[aresta];
  // printf("entrou PhElem::set_border_tipo: border_num=%d\n",en);
  border[en].tipo=t;
  // Aqui termina a marcacao do tipo da borda
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_border_num(const int & aresta, const int & num)
{
  border_num[aresta]=num;
};
// ****************************************************************************
template<int NumVariaveis>
int PhElem<NumVariaveis>::get_trace_border_map(const int & aresta)
{
  int qmax=ptr_stdel[0]->qborder_val();
  int temp = gbtrbmap+qmax*aresta;
  return temp;
}
;
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_type(const int & num)
{
  type = num;
};
// ****************************************************************************
template<int NumVariaveis>
int PhElem<NumVariaveis>::type_val()
{
  return type;
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::print_numeracao(FILE * fout,const int & k)
{
  fprintf(fout,"Nos nas arestas da variavel %d\n",k);
  for(int i=0; i<numb[k];i++){
    fprintf(fout,"print_numeracao nl=%4d ng=%4d\n",i,gbnmap[k][i]);
  }
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::print_modes(FILE * fout,const int & ia) const
{
  int _nv=ptr_stdel[ia]->nv_val(); // number of vertices
  fprintf(fout,"modes da variavel %d = %d",ia,_nv);
  for(int i=0; i<_nv;i++){
    fprintf(fout," %d",gbnmap[ia][i]);
  }
  fprintf(fout,"\n");
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::read_vertices(FILE * finput, const int & _nv) // nao usada
{
  int label;
  for(int i=0;i<_nv;i++)
    fscanf(finput,"%d",&Vert_map[i]);
  fscanf(finput,"%d",&label);
};

// ****************************************************************************
template<int NumVariaveis>
int PhElem<NumVariaveis>::map(const int & var,const int & no_local)
{
  return (gbnmap[var][no_local]);
};

// ****************************************************************************

//template<int NumVariaveis>
//int PhElem<NumVariaveis>::show_NumLocalVars()
//{
//  return(NumLocalVars);
//};

// ***************************************************************************
template<int NumVariaveis>
Stdel * PhElem<NumVariaveis>::show_ptr_stdel(const int & n) const
{
  return ptr_stdel[n];
};

// ****************************************************************************
template<int NumVariaveis>
double PhElem<NumVariaveis>::show_u(const int & var,const int & no)
{
  return u0[var][no];
};

// ****************************************************************************
template<int NumVariaveis>
double PhElem<NumVariaveis>::show_b0(const int & var, const int & no)
{
  return (b0[var][no]);// * sgn[var][no]);
};

// ****************************************************************************
template<int NumVariaveis>
double PhElem<NumVariaveis>::show_bs(const int & var, const int & no)
{
  return (bs[var][no]);// * sgn[var][no]);
};


// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::P_eval_phys(const int & ivar,double f0[])
{
  int nn=numn[ivar];
  int i;
  double utemp[nn];
  for(i=0;i<nn;++i){
    utemp[i]=u0[ivar][i];
  }
  ptr_stdel[ivar]->evalGQ(f0,utemp);
};


// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::P_print(const int & ivar, FILE *file)
{
  int nn=numn[ivar];
  int i;
  double utemp[nn];
  for(i=0;i<nn;i++){
    utemp[i]=u0[ivar][i];
  }
  ptr_stdel[ivar]->printtofile(file,utemp,ptvert,Vert_map);
};


// ****************************************************************************
// Processar dados dos elementos apontando para                               *
// a rotina de cada tipo de elemento                                          *
// ****************************************************************************
/* //se tornou obsoleta
 template<int NumVariaveis>
 void  PhElem::Processar_dados(const int & nel,int& NG,int& NL,
 std::vector<EDGE>& border,
 int Ng[], //inicio da numeracao da borda
 int & NF, int Face[], int Fng[], //inicio da numeracao da face
 int f_mask[])
 {
 int gbnmap_temp[MAXNN]; // mapping from local to global modes (or modes)
 int sgn_temp[MAXNN];
 int sinal_temp[12];
 cout << "PhElem::Processar_dados (velha)\n";
 ptr_stdel[0]->Processar_geometria(nel,ptvert,numv,Vert_map,
 gbnmap_temp,sgn_temp, // candidato a sair
 sinal_temp,
 NG,NL,border,
 Ng, // candidato a sair
 NF,Face,
 Fng,  // candidato a sair aguardar ver se necessario em Tetrahedral
 f_mask);

 //  for(int i=0;i<ptr_stdel[0]->nn_val(); i++){
 //     gbnmap[0][i]=gbnmap_temp[i]; // gbnmap do modo
 //     sgn[0][i]=sgn_temp[i]; // sinal do modo
 //   }

 for(int i=0;i<ptr_stdel[0]->ne_val();i++) sinal[i]=sinal_temp[i]; // sinal da aresta
 //printf("PhElem::Ler PONTO 1c\n");
 };
 */

// **************************************************************************
// Candidato a substituir a funcao acima
//
template<int NumVariaveis>
void PhElem<NumVariaveis>::Processar_dados(int& NL,
                                           std::vector<ARESTA>&aresta_vec,
                                           int & NF,
                                           std::vector<FACE> & face_vec)
{
  int na,nb;
  //cout << "PhElem<NumVariaveis>::Processar_dados (nova)\n";
  // recuperar os numeros de vertices, arestas e faces do elemento padrao
  ndim = ptr_stdel[0]->ndim_val();
  numv=ptr_stdel[0]->nv_val();
  nume=ptr_stdel[0]->ne_val();
  numf=ptr_stdel[0]->nf_val();
  numborders =ptr_stdel[0]->nborder_val();

  // Processar as arestas
  if(ndim!=1) {
    for(int i = 0; i<nume;++i){
      int n0=Vert_map[ptr_stdel[0]->aresta_lvert(i,0)];
      int n1=Vert_map[ptr_stdel[0]->aresta_lvert(i,1)];
      if(n0<n1){ // vertices em ordem crescente
        sinal[i]=1;
        na=n0;
        nb=n1;
      }
      else { // vertices em ordem decrescente
        sinal[i]=-1;
        na=n1;
        nb=n0;
      }
      aresta_map[i] = aresta_gbnum(na,nb,NL,aresta_vec);
      // cout << "aresta "<< aresta_map[i] << " ("<< na<< ","<< nb<< ")"<< endl;
    }
  }

  // Processar as faces
  if(ndim==3) {
    for(int i=0;i<numf;++i){
      int nvf = ptr_stdel[0]->show_nvf(i);
      int var[nvf];
      for(int j=0;j<nvf;++j){
        var[j]= Vert_map[ptr_stdel[0]->face_lvert(i,j)];
      }
      int n = face_gbnum(nvf,var,NF,face_vec);
      face_map[i] = n;
      if(face_vec[n]._tipo == 0){
        face_vec[n]._tipo=ptr_stdel[0]->show_face_tipo(i);
      }
    }
  }
  //printf("PhElem<NumVariaveis>::Ler PONTO 1c\n");
};

// ***************************************************************************
// Salva o numero de vertices e o o mapa de Vertices Vert_map
// ***************************************************************************
// 19/07/2014
/*
 template<int NumVariaveis>
 void PhElem::set_Vert_map(const int & n_in,int ver_temp[])
 {
 numv=n_in;

 if(type==4){ // tetrahedral
 int seq_orig[4] = {0,1,2,3};
 ordenar(4,ver_temp,seq_orig);
 face_mask = new int[4];
 int facemap[4][3] = // sequencia das faces do tetraedro
 { 0,1,2,
 0,1,3,
 0,2,3,
 1,2,3
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
 }
 }
 */

template<int NumVariaveis>
void PhElem<NumVariaveis>::set_Vert_map(const int & n_in,int ver_temp[])
{
  numv=n_in;
  for(int i=0;i<n_in;++i)Vert_map[i]=ver_temp[i];

}
// ***************************************************************************
// Salva os mapas (local para global) e os sinais
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_gbnmap(const int & ia,const int gbnmap_temp[],
                                            const int sgn_temp[])
{
  for(int i=0;i<numn[ia]; i++){
    gbnmap[ia][i]=gbnmap_temp[i];
    sgn[ia][i]=sgn_temp[i];
  }
};

// *********************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::check_connectivity(FILE * fout,const int & ia)
{
  int i,p,q,r;
  for(i=0; i<numb[ia]; ++i){
    ptr_stdel[ia]->show_ind(i,p,q,r);
    fprintf(fout,"i = %3d gbnmap = %3d p = %3d q = %3d r = %3d\n", i, gbnmap[ia][i],p,q,r);
  }
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::check_gradiente(FILE * fileout,const int & ia)
{
  int NGQP=ptr_stdel[ia]->NGQP_val();
  double *ptr_grad[ndim];
  double grad[ndim][NGQP];
  double sn[NGQP];
  for(int i=0;i<ndim;i++){
    ptr_grad[i]=grad[i];
  }
  ptr_stdel[ia]->evalGQ(sn,u0[ia]);
  ptr_stdel[ia]->Gradiente(fileout,ptr_grad,sn,ptvert,Vert_map);

  //ptr_stdel[ia]->Gradiente(fileout,ptr_grad,funcao,ptvert,Vert_map);
  // for(int m=0;m<NGQP;m++)
  //fprintf(fileout,"check_gradiente %11.4e\n",grad[0][m]);
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::assign_gbnmap(const int & ia,const int & i, const int & val)
{
  gbnmap[ia][i]=val;
  //  sgn[ia][i]=1;
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_gbnmap(int & count)
{
  for(int k=0;k<NumVariaveis;++k){
    inicia_gbnmap(k,count);
  }
};
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_gbnmap(const int & ivar,int & count)
{ // Nao impoe continuidade entre os elementos vizinhos
  int NN=numn[ivar];
  for(int i=0;i<NN;++i){
    gbnmap[ivar][i]=count++;
    sgn[ivar][i]=1;
  }
};
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_gbtrbmap(int & count)
{
  // qmax= number of quadrature points on each border
  const int N = numborders;
  const int qmax= ptr_stdel[0]->qborder_val();
  //a ser substituido por numborders(number of borders)
  stgbtrbmapM = new int [N];
  stgbtrbmapP = new int [N];
  gbtrbmap = count;
  count += (N*qmax);
  // printf(" gbtrbmap %d\n", gbtrbmap);
};
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_stgbtrbmap(const int & b,const int & vmapM, const int & vmapP)
{
  stgbtrbmapM[b] = vmapM;
  stgbtrbmapP[b] = vmapP;
};

// ******************************************************************************/
// Pressupoem que o JV (vetor dos Jacobianos nos pontos de Gauss ja foi calculado
// ******************************************************************************/
template<int NumVariaveis>
void PhElem<NumVariaveis>::projetar_C0(FILE *file,
                                       double (*func)(double,double,double),
                                       const int & ivar)
{
  //printf("\nComeco de PhElem::projetar_C0 variavel %d\n",ivar);
  //int i,j,ii,jj,k;
  int _nb=/*numb[ivar];*/ ptr_stdel[ivar]->nb_val();
  int _nn=/*numn[ivar];*/ ptr_stdel[ivar]->nn_val();
  int _ni=_nn-_nb;
  //cout << "PhElem::projetar_C0 ni = "<< ni << endl;
  double aux=0.0;
  int nmap[_nn],bflag[_nn],sgn[_nn];
  double Xl[_nn];// coeficientes dos modos locais; vetor a ser calculado
  //cout << "\n Projetar C0\n";
  // cout << "nb "<< nb << " nn "<< nn <<" nborder "<< nborder<< "\n";
  for(int i=0;i<_nn;i++){
    nmap[i]=i;
    sgn[i]=1;
    bflag[i]=1;
    Xl[i]=0.0;
  }

  for(int j=0;j<numborders;j++){
    cout << "Chamar Dirichlet (ptr_stdel) para border = "<< j << endl;
    ptr_stdel[ivar]->Dirichlet(j,ptvert,Vert_map,nmap,sgn,bflag,Xl,func);
    // cout << "saindo de Dirichlet para a face "<< j << endl;
    //for(int k=0;k<nn;++k) cout << "Xl["<< k << "%d] = " <<Xl[k] << endl;
  }
  cout << "Terminou loop sobre as bordas (ptr_stdel)\nRetornou a PhElem::projetar\n";

  //ni=0;
  if(_ni>0){
    //cout << "calcular o produto interno de f por phi dos "<< ni << " modos internos"<< endl;

#ifdef _NEWMAT
    // newmat
    NEWMAT::Matrix Mi(_ni,_ni);
    NEWMAT::ColumnVector B(_ni);
#endif

    //fprintf(file,"\n\n%5d\n",ni*ni);
    //printf("Matriz Mi\n");
    for(int i=_nb;i<_nn;++i){
      int ii=i-_nb;
      Mi.element(ii,ii)=ptr_stdel[ivar]->mass(i,i,JV);
      // printf("%3d %3d %g\n",ii,ii,Mi.element(ii,ii));
      for(int j=i+1;j<_nn;j++){
        int jj=j-_nb;
        aux=ptr_stdel[ivar]->mass(i,j,JV);
        Mi.element(ii,jj)=aux;
        Mi.element(jj,ii)=aux;
        //printf("%3d %3d %g\n%3d %3d %g\n",ii,jj,aux,jj,ii,aux);
      }
    }

    int _NGQP = ptr_stdel[ivar]->NGQP_val();
    double phi[_NGQP];
    double f[_NGQP];
    // inicializar f;
    ptr_stdel[ivar]->computeFuncGQ(f,ptvert,Vert_map,func);
    for(int j=0;j<_nb;j++){
      ptr_stdel[ivar]->eval_Phi(j,phi);
      double aux = Xl[j];
      for(int k=0;k<_NGQP;k++) {f[k] -= (aux*phi[k]);}
    }
    double b[_nn];
    ptr_stdel[ivar]->vector_of_integral_of_f_Phi_dv(b,f,JV);
    //fprintf(file,"%5d\n",ni);
    for(int i=0;i<_ni;i++){
      B.element(i)=b[i+_nb];
      //fprintf(file,"%g ",B[i]);
    }
    //fprintf(file,"\n");

#ifdef _NEWMAT
    NEWMAT::ColumnVector Y = Mi.i() * B; // B=Y ; // newmat
#endif

    for(int i=_nb;i<_nn;++i)u0[ivar][i]=Y.element(i-_nb);
  }
  // fim de if(ni>0)

  for(int i=0;i<_nb;++i) u0[ivar][i]=Xl[i]; // no contorno

  if(file != NULL)
    ptr_stdel[ivar]->printtofile(file,u0[ivar],func,ptvert,Vert_map);

  cout << "Terminou  PhElem::projetar_C0\n\n";

};

// ***********************************************************************
// calcula os coeficientes a partir dos valores nos pontos de Gauss      *
// ***********************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::transformacao_direta(double f[],const int & ivar)
{
  int i,j;
  int nn=numn[ivar];//ptr_stdel[ivar]->nn_val();
  double aux;
  //cout << "\n Projetar C0\n";

#ifdef _NEWMAT
  // newmat
  NEWMAT::Matrix Mi(nn,nn);
  NEWMAT::ColumnVector B(nn);
#endif

  //fprintf(file,"\n\n%5d\n",ni*ni);
  for(i=0;i<nn;i++){
    Mi.element(i,i)=ptr_stdel[ivar]->mass(i,i,JV);
    //fprintf(file,"%3d %3d %g\n",ii,ii,Mi[ii][ii]);
    for(j=i+1;j<nn;j++){
      aux=ptr_stdel[ivar]->mass(i,j,JV);
      Mi.element(i,j)=aux;
      Mi.element(j,i)=aux;
      //fprintf(file,"%3d %3d %g\n%3d %3d %g\n",ii,jj,aux,jj,ii,aux);
    }
  }

  double b[nn];
  ptr_stdel[ivar]->vector_of_integral_of_f_Phi_dv(b,f,JV);
  for(i=0;i<nn;i++){
    B.element(i)=b[i];
  }
#ifdef _NEWMAT
  NEWMAT::ColumnVector Y = Mi.i() * B;  B=Y ; // newmat
#endif

  for(i=0;i<nn;i++)u0[ivar][i]=B.element(i);
};

// ****************************************************************************
// template<int NumVariaveis>
// void PhElem<NumVariaveis>::BoundCond(const int & face,int bflag[],double X[],
// 		     double (*func)(double,double,double),
// 		     const int & nf, const int & bndtype)
// {
//   // ***********************************************************************
//   // bndtype = 0 Dirichlet; marca bflag = 0 (variavel conhecida)
//   // bndtype = 1 Especifica forca; marca bflag = 1 (variavel desconhecida)
//   // ***********************************************************************
//   if(bndtype==0)
//     ptr_stdel[nf]->Dirichlet(face,ptvert,gbnmap[nf],sgn[nf],bflag,X,func,nf,NumLocalVars);
//   if(bndtype==1){
//     printf("BNDTYPE = 1\n");
//     ptr_stdel[nf]->BoundForce(face,ptvert,gbnmap[nf],bs[nf],func,nf,NumLocalVars);
//     }
// }


//  // ************************************************************************
//  // ************************************************************************
//  template<int NumVariaveis>
//  void PhElem<NumVariaveis>::AlocarMatrizesK(Matrix * pt1,Matrix * pt2,Matrix *pt3)
//  {
//    ptKcomp=pt1;
//    ptKi_inv=pt2;
//    ptKcKi_inv=pt3;
//  };

//    // **********************************************************************
//    template<int NumVariaveis>
//    void PhElem<NumVariaveis>::AlocarMatrizesM(Matrix * pt1,Matrix * pt2,Matrix *pt3)
//    {
//      ptMb_comp=pt1;
//      ptMi_inv=pt2;
//      ptMcMi_inv=pt3;
//    };
//    // **********************************************************************
//    template<int NumVariaveis>
//    void PhElem<NumVariaveis>::make_MassMatrices()
//    {
//      int i,j,ii,jj,im,jm;
//      int nb=ptr_stdel->nb_val();
//      int Nb=nb*NumLocalVars;
//      int nn=ptr_stdel->nn_val();
//      int Nn=nn*NumLocalVars;
//      int Ni=Nn-Nb;
//      int ni=nn-nb;
//
//      // Alocacao dinamica das matrizes
//      ptMb_comp   = new Matrix(0,Nb-1,0,Nb-1);
//      ptMi_inv  = new Matrix(0,Ni-1,0,Ni-1);
//      ptMcMi_inv= new Matrix(0,Nb-1,0,Ni-1);
//
//      for(int i=0;i<Nb;i++){
//        for(int j=0;j<Nb;j++)
//          (*ptMb_comp)(i,j)=0.0;
//        for(int j=0;j<Ni;j++)
//          (*ptMcMi_inv)(i,j)=0.0;
//      }
//      for(int i=0;i<Ni;i++)
//        for(int j=0;j<Ni;j++)
//          (*ptMi_inv)(i,j)=0.0;
//
//      NEWMAT::Matrix mb(0,nb-1,0,nb-1);
//      NEWMAT::Matrix mb_c(0,nb-1,0,nb-1);
//      NEWMAT::Matrix mi(0,ni-1,0,ni-1);
//      NEWMAT::Matrix mc(0,nb-1,0,ni-1);
//      NEWMAT::Matrix mi_inv(0,ni-1,0,ni-1);
//      NEWMAT::Matrix mcmi_inv(0,nb-1,0,ni-1);
//
//      // Matriz Mb
//      for(i=0;i<nb;i++){
//        for(j=0;j<nb;j++) mb.element(i,j)=ptr_stdel->mass(i,j,JV);
//        // Matrix Mc
//        for(j=nb;j<nn;j++){
//          jj=j-nb;
//          mc.element(i,jj)=ptr_stdel->mass(i,j,JV);
//        }
//      }
//      // Matrix Mi
//      for(i=nb;i<nn;i++){
//        ii=i-nb;
//        mi.element(ii,ii)=ptr_stdel->mass(i,i,JV);
//        for(j=i+1;j<nn;j++){
//          jj=j-nb;
//          double aux=ptr_stdel->mass(i,j,JV);
//          mi.element(ii,jj)=aux;
//          mi.element(jj,ii)=aux;
//        }
//      }
//
//      mi_inv=mi.Inverse();
//      mcmi_inv = mc * mi_inv;
//      mb_c = mb - (mcmi_inv * mc.Transpose());
//
//      for(int ivar=0;ivar<NumLocalVars;ivar++){
//        for(i=0;i<nb;i++){
//          im=i*NumLocalVars+ivar;
//          for(j=0;j<nb;j++)
//    	(*ptMb_comp)[im][j*NumLocalVars+ivar]=mb_c.element(i,j);
//          for(j=nb;j<nn;j++){
//    	jj=j-nb;
//    	jm=jj*NumLocalVars+ivar;
//    	(*ptMcMi_inv)[im][jm]=mcmi_inv.element(i,jj);
//          }
//        }
//        //printf("Making local matrices 3\n");
//        for(i=nb;i<nn;i++){
//          ii=i-nb;
//          im=ii*NumLocalVars+ivar;
//          for(j=nb;j<nn;j++){
//    	jj=j-nb;
//    	jm=jj*NumLocalVars+ivar;
//    	(*ptMi_inv)[im][jm]=mi_inv.element(ii,jj);
//          }
//        }
//      }
//    }


// ****************************************************************************
// PhElem<NumVariaveis>::Vector(int i) creates the vector for the boundary modes performing
// the Schur condensation Vector = fb-McMi_inv*fi and applying the sgn[i]
// ****************************************************************************
// template<int NumVariaveis>
// double PhElem<NumVariaveis>::Vector(const int & ia,const int & i)
// {
//     int nb=ptr_stdel[ia]->nb_val();
//     int ni=(ptr_stdel[ia]->nn_val()-nb);
//     double sum0=b0[ia][i];
//     for(int j=0; j<ni; j++)
//       sum0-=(*ptMcMi_inv)(i,j)*b0[ia][j+nb];
//     return (sum0*sgn[ia][i]);
// };

//  // ***********************************************************************
// template<int NumVariaveis>
//  void PhElem<NumVariaveis>::Vector2(const int & i,double sum[])
//  {
//    // *********************************************************
//    // i = numero local do no; cada no tem NumLocalVars variaveis *
//    // *********************************************************
//    //printf("entrou em vector2\n");
//    int nb=ptr_stdel->nb_val();// numero de nos de contorno
//    int ni=ptr_stdel->nn_val()-nb;// numero de nos internos
//    int Ni=ni*NumLocalVars;
//    int Nb=nb*NumLocalVars;
//    int ia;
//    for(ia=0;ia<NumLocalVars;ia++){
//      int ii = NumLocalVars*i + ia;
//      sum[ia]=b0[ii];
//      for(int j=0; j<Ni; j++){
//        sum[ia]-=(*ptMcMi_inv)(ii,j)*b0[j+Nb];
//      }
//    }
//    for(ia=0;ia<NumLocalVars;ia++)
//      sum[ia]*=sgn[i];
//  };

//  // ************************************************************************
//  double PhElem<NumVariaveis>::MassMatrix_value(int i,int j)// <== i e j sao modos locais=<
//  {
//    //printf("entrou em MassMatrix_value\n");
//    return((*ptMb_comp)(i,j)*sgn[(i/NumLocalVars)]*sgn[(j/NumLocalVars)])*rho;
//  };
//

// ****************************************************************************
// int PhElem<NumVariaveis>::map(const int & no_local)
// {
//   return (gbnmap[no_local]);
// };

//  // ************************************************************************
//  void PhElem<NumVariaveis>::P_eval_u(const double X[])
//  {
//    int nb=ptr_stdel->nb_val();
//    int Nb=nb*NumLocalVars;
//    int Nn=ptr_stdel->nn_val()*NumLocalVars;
//    int ii,gbi;
//
//    for(int ia=0;ia<NumLocalVars;ia++){
//      for(int i=0;i<nb;i++){
//        ii=i*NumLocalVars+ia;
//        gbi=map(i,ia);
//        u0[ii]=X[gbi]*sgn[i];
//      }
//    }
//    for(int i=Nb;i<Nn;i++){
//      u0[i]=0.0;
//      for(int j=0;j<Nb;j++){
//        u0[i]-=(*ptMcMi_inv)(j,i-Nb)*u0[j];
//      }
//      for(int j=Nb;j<Nn;j++){
//        u0[i]+=(*ptMi_inv)(i-Nb,j-Nb)*b0[j]/rho;
//      }
//    }
//  };
// ******************************************************************
/* void  PhElem::P_eval_u(const double X[],const int & ia)
 {
 int nb=ptr_stdel[ia]->nb_val();
 int nn=ptr_stdel[ia]->nn_val();
 int gbi;
 for(int i=0;i<nb;i++){
 gbi=map(i,ia);
 u0[ia][i]=X[gbi]*sgn[ia][i];
 }
 for(int i=nb;i<nn;i++){
 u0[ia][i]=0.0;
 for(int j=0;j<nb;j++){
 u0[ia][i]-=(*ptMcMi_inv)(j,i-nb)*u0[ia][j];
 }
 for(int j=nb;j<nn;j++){
 u0[ia][i]+=(*ptMi_inv)(i-nb,j-nb)*b0[ia][j]/rho;
 }
 }
 };
 */
// // ******************************************************************
// void PhElem<NumVariaveis>::P_eval_u_p1(const double X[],const int & ia)
// {
//   // Valida para a Equacao_p1 do problema ImiscFlow2F
//   int nb=ptr_stdel[ia]->nb_val();
//   int nn=ptr_stdel[ia]->nn_val();
//   int ii,jj,gbi;
//   for(int i=0;i<nb;i++){
//     gbi=map(i,ia);
//     u0[ia][i]=X[gbi]*sgn[ia][i];
//   }
//   for(int i=nb;i<nn;i++){
//     u0[ia][i]=0.0;
//     for(int j=0;j<nb;j++){
//       u0[ia][i]-=(*ptKcKi_inv)(j,i-nb)*u0[ia][j];
//     }
//     for(int j=nb;j<nn;j++){
//       u0[ia][i]+=(*ptKi_inv)(i-nb,j-nb)*b0[ia][j];
//     }
//   }
// };


// ****************************************************************************
//int PhElem<NumVariaveis>::show_sgn(const int & ia, const int & i)
//{
// return (sgn[ia][i]);
// };

//  double PhElem<NumVariaveis>::show_u(const int & ia,const int & i)
//  {
//    return u0[ia][i];
//  };

// ****************************************************************************
// Calcula os valores da funcao de interpolacao nos vertices
// ****************************************************************************
//  void PhElem<NumVariaveis>::P_eval_vert(const double X[],double f_vert[])
//  {
//    int i,m1;
//    int nb=ptr_stdel->nb_val();
//    int nn=ptr_stdel->nn_val();
//    int ni=nn-nb;
//    int p, q;
//    for(i=0;i<nb;i++)u0[i]=X[gbnmap[i]]*sgn[i];
//    for(i=nb;i<nn;i++){
//      u0[i]=0.0;
//      for(int j=0;j<nb;j++){
//        u0[i]-=(*ptMcMi_inv)(j,i-nb)*u0[j];
//      }
//      for(int j=nb;j<nn;j++){
//        u0[i]+=(*ptMi_inv)(i-nb,j-nb)*b0[j]/rho;
//      }
//    }
//    ptr_stdel->computeVertice(f_vert,u0,ptvert,gbnmap);
//  };

// ***************************************************************************
// void PhElem<NumVariaveis>::P_eval_print(const double X[],const int & numf, FILE *file,
// 			double (*func)(double,double,double))
// {
//   int nn=ptr_stdel[numf]->nn_val();
//   int i;
//   double utemp[nn];
//   P_eval_u(X,numf);
//   for(i=0;i<nn;i++){
//     utemp[i]=u0[numf][i];
//     //printf("P_eval_print var=%d u[%d]=%g\n",numf,i,utemp[i]);
//   }
//   ptr_stdel[numf]->printtofile(file,utemp,func,ptvert,Vert_map);
// };


// ***************************************************************************
//  void PhElem<NumVariaveis>::print_matrices(FILE * fout)
//  {
//    int i,j;
//    int nb=ptr_stdel->nb_val();
//    int nn=ptr_stdel->nn_val();
//
//    // fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mass\n");
//    // for(i=0;i<Nb;i++){
//    //   //fprintf(fout,"linha %d\n",i);
//    //   for(j=0;j<Nb;j++){
//    //     fprintf(fout,"%11.4e ",ptr_stdel->mass(i,j)*rho);
//    //   }
//    //   fprintf(fout,"\n");
//    // }
//
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mb_comp\n");
//    for(i=0;i<nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      for(j=0;j<nb;j++){
//        fprintf(fout,"%11.4e ",(*ptMb_comp)(i,j)*rho);
//      }
//      fprintf(fout,"\n");
//    }
//
//    //  fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mc_transposta\n");
//    //   for(i=0;i<nn-nb;i++){
//    //     //fprintf(fout,"linha %d\n",i);
//    //     for(j=0;j<nb;j++){
//    //       fprintf(fout,"%11.4e ",Mc_T[i][j]);
//    //     }
//    //     fprintf(fout,"\n");
//    //   }
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mi_inv\n");
//    for(i=0;i<nn-nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      for(j=0;j<nn-nb;j++){
//        fprintf(fout,"%11.4e ",(*ptMi_inv)(i,j)/rho);
//      }
//      fprintf(fout,"\n");
//    }
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz McMi_inv\n");
//    for(i=0;i<nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      for(j=0;j<nn-nb;j++){
//        fprintf(fout,"%11.4e ",(*ptMcMi_inv)(i,j));
//      }
//      fprintf(fout,"\n");
//    }
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices vetor b\n");
//    for(i=0;i<nn;i++){
//      //fprintf(fout,"linha %d\n",i);
//      fprintf(fout,"%11.4e ",b0[i]);
//    }
//    fprintf(fout,"\n");
//
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices gbnmap e sgn  \n");
//    for(i=0;i<nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      fprintf(fout,"%4d %4d %4d\n",i,gbnmap[i], sgn[i]);
//    }
//    fprintf(fout,"\n");
//  };

// *********************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::compute_JV(const int & ia)
{
    cout << " PhElem<NumVariaveis>::compute_JV(const int & ia): ia =  " << ia << endl;
    cout << " PhElem<NumVariaveis>::compute_JV: NumVariaveis =  " << NumVariaveis << endl;
  if(ia < NumVariaveis) {
      printf("PhElem<NumVariaveis>::compute_JV Calculo do Jacobiano: variavel %d\n",ia);
    int q0=ptr_stdel[ia]->NGQP_val();
    
    // Armazenamento dinamico da memoria para JV;
    // Opcao 1: Uma matriz com 3 indices
    //   JV = new double ** [q0];
    //   for(int i=0;i<q0;i++){
    //     JV[i] = new double * [q1];
    //     for(int j=0;j<q1;j++){
    //       JV[i][j] = new double  [q2];
    //     }
    //   }

    //JV[ia] = new double [q0*q1*q2]; // opcao 2: um unico vetor
    JV = new double [q0];

    printf("ANTES: Calculo do Jacobiano da var %d  dimensao=%d\n",ia,q0);

    ptr_stdel[ia]->Jacobian(ptvert,Vert_map,JV);
  }
  else cout << " elemento nao contem a variavel " << ia << endl;
  //printf("DEPOIS: Calculo do Jacobiano dimensao=%d\n",q0*q1*q2);
};


// ******************************************************************************
// ******************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_part_num(const int & num)
{
  if(num < 0) { //calcula
    // ****************************************************
    // Elemento pertence a particao de menor numero
    // dentre as particoes de seus vertices
    // ****************************************************
    int temp;
    part_num=ptvert[Vert_map[0]].part_num;
    for(int i=1;i<numv;++i){
      temp=ptvert[Vert_map[i]].part_num;
      if(part_num > temp) part_num=temp;
    }
  }
  else // part_num eh dado na chamada do metodo
    part_num=num;
};

// *****************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::teste_gradiente()
{
  const int nn=numn[0];//ptr_stdel[0]->nn_val();
  const int NGQP=ptr_stdel[0]->NGQP_val();
  // const int ns=ptr_stdel[sat]->nn_val();
  // const int np=ptr_stdel[pres]->nn_val();
  double ** der = new  double * [ndim];
  double ** grad = new  double * [ndim];
  for(int k=0;k<ndim;k++) {
    der[k] = new double[NGQP];
    grad[k] = new double[NGQP];
  }
  for(int i=0;i<2;i++){
    for(int j=0; j < nn; j++) {// j=modo
      // calcular os valores de Phi_j nos pontos de quadratura
      double phi[NGQP];
      ptr_stdel[i]->eval_Phi(j,phi);
      // Calcular as derivadas de Phi_j
      ptr_stdel[i]->Gradiente(grad,phi,ptvert,Vert_map);
      ptr_stdel[i]->eval_GradPhi(ptvert,Vert_map,j,der);
      for(int q=0;q<NGQP;q++){
        for(int ndir=0;ndir<ndim;ndir++){
          double aux= grad[ndir][q] - der[ndir][q];
          if(abs(aux) > 1e-14 ){
            printf("var = %2d modo = %2d ponto = %2d dir = %2d  grad = %16.9e der = %16.9e diferenca = %16.9e\n",i,j,q,ndir,grad[ndir][q],der[ndir][q], aux);
          }
        }
      }
    }
  }
  for(int k=0;k<ndim;k++) {
    delete [] der[k]; der[k]=nullptr;
  }
  delete [] der; der=nullptr;
};
/*
 template<int NumVariaveis>
 void PhElem::construir_PhiArray()
 {
 int NGQP = ptr_stdel[0]->NGQP_val();
 int nn =  ptr_stdel[0]->nn_val();

 double phi[NGQP];
 PhiArray = new double [nn*NGQP];

 int count = 0;
 for(int j=0; j < nn; j++) {// j=modo
 // calcular os valores de Phi_j nos pontos de quadratura
 ptr_stdel[0]->eval_Phi(j,phi);
 for(int i=0;i<NGQP;i++) PhiArray[count++]=phi[i];
 }
 };
 template<int NumVariaveis>
 const double * PhElem::eval_Phi(const int i)
 {
 int NGQP = ptr_stdel[0]->NGQP_val();
 return (PhiArray + (i*NGQP));
 };
 */

template<int NumVariaveis>
void PhElem<NumVariaveis>::vetor_superficie(const int & num_local,double & area, double normal[3])
{
  ptr_stdel[0]->superficie_externa(ptvert,Vert_map,num_local,area,normal);
};

//   // ****************************************************************************
//  // ****************************************************************************
// template<int NumVariaveis>
//  void PhElem<NumVariaveis>::VolumeIntegrals_IG_Epetra(Fluids fls,
//  		       Epetra_FECrsMatrix & A,
//  		       Epetra_FEVector & RHS)
//  {
//   cout << "PhElem<NumVariaveis>::VolumeIntegrals_IGa \n";
//  #define indice(n0,var,i) ((n0*var)+i)
//  #define npos(ntot,i,j) ((ntot*i)+j)
//
//
//    int qmax=ptr_stdel[sat]->Q_val(0);
//    int NGQP=ptr_stdel[sat]->Q_val(0)*ptr_stdel[sat]->Q_val(1);
//
//    int ns,np;
//    int h,i,j,k,m,ivar;
//    double aux,aaux;
//    double res0[NGQP],res1[NGQP],res2[NGQP],res3[NGQP];
//    double JW[NGQP];
//    double sn[NGQP], pw[NGQP],pc[NGQP];
//    double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
//    double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
//    double d_pc[NGQP],d2_pc[NGQP];
//    double phi_r[NGQP],phi_l[NGQP];
//    double grad_sn[3][NGQP],grad_pw[3][NGQP],grad_pc[3][NGQP];
//    double *ptr_gsn[3],*ptr_gpw[3],*ptr_gpc[3];
//    double mun=fls.show_mun();
//    double muw=fls.show_muw();
//    for(k=0;k<3;k++){
//      ptr_gsn[k]=grad_sn[k];
//      ptr_gpw[k]=grad_pw[k];
//      ptr_gpc[k]=grad_pc[k];
//    }
//
//    int ndim=ptr_stdel[sat]->ndim_val();
//    int ne  =ptr_stdel[sat]->ne_val();
//
//    // calcula os pesos de Gauss multiplicados pelos jacobianos
//    if(ptr_stdel[sat]->w_vec(JW) != NGQP){
//      printf("Erro de dimensionamento de JW\n");
//      exit(0);
//    }
//    for(m=0;m<NGQP;m++)JW[m]*=JV[m];
//
//    // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
//    ptr_stdel[sat]->evalGQ(sn,u0[sat]);
//    ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
//    ptr_stdel[pres]->evalGQ(pw,u0[pres]);
//    ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
//
//    // Calcula pc e seu gradiente nos pontos de Gauss
//    for(m=0;m<NGQP;m++)pc[m]=fls.pressao_capilar(sn[m]);
//    ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
//
//    for(m=0;m<NGQP;m++){
//      // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
//      for(k=0;k<ndim;k++){// caso especial de tensor diagonal
//        grad_sn[k][m]*=perm[k];
//        grad_pw[k][m]*=perm[k];
//        grad_pc[k][m]*=perm[k];
//      }
//      // Calculo dos escalares
//      aux=sn[m];
//      lambdan[m]=fls.Krn(aux)/mun;
//      lambdaw[m]=fls.Krw(aux)/muw;
//      lambdat[m]=lambdan[m]+lambdaw[m];
//      d_lambdan[m]=fls.dkrn(aux)/mun;
//      d_lambdaw[m]=fls.dkrw(aux)/muw;
//      d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
//      d_pc[m]=fls.dpc(aux);
//      d2_pc[m]=fls.d2pc(aux);
//    }
//
//    // salva os tracos de K * gradientes
//    /*
//    for(h=0;h<ne;h++){
//      ptr_stdel[sat]fill_trace(h,qmax,sinal[h],sn,Trsn[h]);
//      ptr_stdel[pres]fill_trace(h,qmax,sinal[h],pw,Trpw[h]);
//      for(k=0;k<ndim;k++){
//        ptr_stdel[sat]fill_trace(h,qmax,sinal[h],grad_sn[k],TrKgrad_sn[h][k]);
//        ptr_stdel[sat]fill_trace(h,qmax,sinal[h],grad_pc[k],TrKgrad_pc[h][k]);
//        ptr_stdel[pres]fill_trace(h,qmax,sinal[h],grad_pw[k],TrKgrad_pw[h][k]);
//      }
//    }
//    */
//    calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
//    np=ptr_stdel[pres]->nn_val();
//    ns=ptr_stdel[sat ]->nn_val();
//    for(int rp=0;rp<np;rp++){ // PE_VI
//      // *****************************
//      // calculo do vetor
//      // *****************************
//
//      // integral 2 com sinal trocado
//      prodgg(NGQP,grad_pc[0],grad_pc[1],
//  	   GradPhi[pres][rp][0],GradPhi[pres][rp][1],res0);
//      prodvv(NGQP,lambdan,res0,res0);
//      integral(NGQP,JW,res0,aux);
//      //B[gbnmap[pres][rp]] -= aux; // primeiro valor de B
//      B[rp] -= aux;
//
//      // Fim do vetor em PE_VI
//
//      // *****************************
//      // Matriz
//      // *****************************
//      for(int lp=0;lp<np;lp++){ // PE_VI_PP
//        prodev(NGQP,perm[0],GradPhi[pres][lp][0],res0);
//        prodev(NGQP,perm[1],GradPhi[pres][lp][1],res1);
//        prodgg(NGQP,res0,res1,GradPhi[pres][rp][0],GradPhi[pres][rp][1],res0);
//        prodvv(NGQP,lambdat,res0,res0);
//        integral(NGQP,JW,res0,aux);
//        //Ti[count]=gbnmap[pres][rp];
//        //Tj[count]=gbnmap[pres][lp];
//        mx[ ( (rp * np) + lp) ] = aux;
//        count++;
//      }
//    }
//    for(int rp=0;rp<np;rp++){
//      indx[rp] = gbnmap[pres][rp];
//    }
//
//    A.InsertGlobalValues(ntot,indx,mx,Epetra_FECrsMatrix::ROW_MAJOR);
//    RHS.SumIntoGlobalValues(ntot,indx,B);
//    delete [] indx;indx=0;
//    delete [] mx; mx=0;
//    delete [] B; B=0;
//    // OK 22/05/2008
//  };



// ***************************************************************************
// Salva as rotinas abaixo para servir de template para casos de
// inicializacao de vetores mais simples que o DG_PhElem
/*
 void PhElem::inicia_vetores()
 {
 for (int var=0; var < NumVariaveis; ++var)
 {
 int Nn = pt_stdel[var]->nn_val();
 b0[var] = new double [Nn];
 u0[var] = new double [Nn];
 //u1 = new double [Nn];
 //u2 = new double [Nn];
 for(int i=0;i<Nn;i++){
 b0[var][i]=0.0;
 u0[var][i]=0.0;
 // u1[i]=0.0;
 //u2[i]=0.0;
 }
 int Nb=pt_stdel[var]->nb_val();
 bs[var] = new double [Nb];
 for(int i=0;i<Nb;i++)bs[var][i]=0.0;
 }
 };
 // *************************************************************
 void PhElem::finaliza_vetores()
 {
 for (int var=0; var < NumVariaveis; ++var)
 {
 delete [] b0[var];  b0[var] = nullptr;
 delete [] u0[var];  u0[var] = nullptr;
 delete [] bs[var];  bs[var] = nullptr;
 }
 delete [] b0; b0 = nullptr;
 delete [] u0; u0 = nullptr;
 delete [] bs; bs = nullptr;
 };

 // *************************************************************
 void PhElem::zera_vetores()
 {
 for (int var=0;var < NumVariaveis; ++var)
 {
 int Nn=pt_stdel[var]->nn_val();
 for(int i=0;i<Nn;i++){
 b0[var][i]=0.0;
 u0[var][i]=0.0;
 //u1[i]=0.0;
 //u2[i]=0.0;
 }
 int Nb=pt_stdel[var]->nb_val();
 for(int i=0;i<Nb;i++)bs[var][i]=0.0;
 }
 };
 void PhElem::zera_bs()
 {
 for (int var=0;var < NumVariaveis; ++var)
 {
 int Nb=pt_stdel[var]->nb_val();
 for(int i=0;i<Nb;i++)bs[var][i]=0.0;
 }
 };
 */

// *****************************************************************************
//#include "GBNMAP.cpp"

// *****************************************************************************
// Mapeamento dos vertices
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_vertices(std::vector< std::vector<int>> gbn)//int ** gbn )
{
  for(int i=0;i<NumVariaveis;++i){
    for(int j=0;j<numv;++j){
      gbnmap[i][j]=gbn[i][Vert_map[j]];
      sgn[i][j]=1;
    }
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_vertices(int gbn[],const int & ivar)
{
  for(int i=0;i<numv;++i){
    gbnmap[ivar][i]=gbn[Vert_map[i]];
    sgn[ivar][i]=1;
  }
};

// *****************************************************************************
// Mapeamento das arestas
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_aresta(const int & aresta,const int & sinal,
                                               const int & ivar,const int & inicio,
                                               int & inc)
{
  int ip;
  int ini = ptr_stdel[ivar]->show_emapi(aresta)+1;
  int fim = ptr_stdel[ivar]->show_emapi(aresta+1)-1;
  int aux=1;
  //cout << "(inicio,fim) = (" << ini << ","<< fim << ")   inc_inicial " << inc << endl;
  for(int i = ini; i < fim; ++i){
    ip=ptr_stdel[ivar]->show_emapv(i);
    gbnmap[ivar][ip]=inicio+inc;
    sgn[ivar][ip]=aux;
    aux*=sinal; // muda o sinal de aux se sinal=-1
    ++inc;
  }
  //cout << "inc_final "<< inc << endl;
};

// *****************************************************************************
// Mapeamento dos modos interiores
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_interior(const int & ivar, int & count)
{
  int ini=numb[ivar];
  int fim=numn[ivar];
  for(int i=ini;i<fim;++i){
    gbnmap[ivar][i]=count++;
    sgn[ivar][i]=1;
  }
};

// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_faces(const int face_vec[],const int & ivar)
{
  int i[3];
  cout << "PhElem<NumVariaveis>::gbnmap_faces\n";
  // Loop sobre as faces
  for(int f=0;f<numf;++f){

    int dir[3];
    int v2=0;
    int sinal[2];

    // achar os vertices da face
    int nvf = ptr_stdel[ivar]->show_nvf(f);
    int Av[nvf], ind[nvf];
    for(int i = 0; i < nvf; ++i){
      ind[i]=ptr_stdel[ivar]->face_lvert(f,i);
      Av[i]=Vert_map[ind[i]];
    }

    if(nvf == 4) {
      quad_ordem(Av,ind,dir,sinal,v2);

      if(v2 == 0) i[dir[2]] = 0; // indice na direcao perpendicular
      else i[dir[2]] = ptr_stdel[ivar]->P_val(dir[2]); // indice na direcao perpendicular
      int P0 = ptr_stdel[ivar]->P_val(dir[0]);
      int P1 = ptr_stdel[ivar]->P_val(dir[1]);

      int inicio = face_vec[face_map[f]];
      int inc = 0;
      int aux0 = 1;
        // numerar modos internos da face
      for(int i0 = 1; i0 < P0; ++i0) {
        i[dir[0]] = i0;
        int aux1 = 1;
        for(int i1 = 1; i1< P1; ++i1) {
          i[dir[1]] = i1;
          int n = ptr_stdel[ivar]->show_ind_mode(i[0],i[1],i[2]);
          gbnmap[ivar][n]=inicio + inc;
          sgn[ivar][n]=aux1*aux0;
          inc++;
          aux1 *= sinal[1];
        }
        aux0 *= sinal[0];
      }
    }
    else if ( nvf == 3) { // faces triangulares
      dir[0]=ptr_stdel[ivar]->show_fd0(f);
      dir[1]=ptr_stdel[ivar]->show_fd1(f);
      dir[2]=3-dir[0]-dir[1];
      v2=ptr_stdel[ivar]->show_fv2(f);

      if(v2 == 0) i[dir[2]] = 0;
      else i[dir[2]] = ptr_stdel[ivar]->P_val(dir[2]);
      int P0 = ptr_stdel[ivar]->P_val(dir[0]);
      int P1 = ptr_stdel[ivar]->P_val(dir[1]);

      int inicio = face_vec[face_map[f]];
      int inc = 0;
      //int aux0 = 1;
        
        // numerar modos internos da face
      for(int i0 = 1; i0 < P0; ++i0) {
        i[dir[0]] = i0;
        //int aux1 = 1;
        for(int i1 = 1; i1< P1 - i0; ++i1) {
          i[dir[1]] = i1;
          int n = ptr_stdel[ivar]->show_ind_mode(i[0],i[1],i[2]);
          gbnmap[ivar][n]=inicio + inc;
          sgn[ivar][n]=1;//aux1*aux0;
          inc++;
          //aux1 *= sinal[1];
        }
        //aux0 *= sinal[0];
      }
    }
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_arestas(const int aresta_vec[],const int & ivar)
{
  int av[2];
  // Loop sobre as arestas
  for(int a=0;a<nume;++a){
    av[0]=ptr_stdel[ivar]->aresta_lvert(a,0);
    av[1]=ptr_stdel[ivar]->aresta_lvert(a,1);
    int sinal = 1;
    if(Vert_map[av[0]] > Vert_map[av[1]]) sinal = -1;

    int ini = ptr_stdel[ivar]->show_emapi(a)+1;
    int fim = ptr_stdel[ivar]->show_emapi(a+1)-1;
    int inicio = aresta_vec[aresta_map[a]];
    int aux=1;
    int inc = 0;
    for(int i = ini; i < fim; ++i){
      int ip=ptr_stdel[ivar]->show_emapv(i);
      gbnmap[ivar][ip]=inicio+inc;
      sgn[ivar][ip]=aux;
      aux*=sinal; // muda o sinal de aux se sinal=-1
      inc++;
    }
  }
};
// ***********************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_vetores()
{
  int _nn;//nb,q0,q1;
  int i,k;//h, j
 cout<< "PhElem<NumVariaveis>::inicia_vetores()\n";
  if(vetores_iniciados == 1) { cout<< "vetores locais de PhElem já iniciados\n"; }

  else {
    vetores_iniciados = 1;
    cout<< "PhElem<NumVariaveis>::inicia_vetores()\n";
    for (k=0;k<NumVariaveis;++k){
      _nn=numn[k];
      u0[k] = new double [_nn];
      usave[k] = new double [_nn];
      for(i=0;i<_nn;++i){
        u0[k][i]=0.0;
      }
    }
    // **************************************
    // Calcular o Jacobiano   *             *
    // **************************************
    // q0= ptr_stdel[0]->Q_val(0);
    // q1= ptr_stdel[0]->Q_val(1);

    //int ndim =  ptr_stdel[0]->ndim_val();
    //int qmax = ptr_stdel[0]->qborder_val();
    //int nborder = ptr_stdel[0]->nborder_val();
    //int NGQP;
    //printf("Calculo do Jacobiano: dimensao=%d\n",NGQP);
    //JV = new double [NGQP]; // opcao 2: um unico vetor
    //ptr_stdel[0]->Jacobian(ptvert,Vert_map,JV);
      cout<< "PhElem<NumVariaveis>::inicia_vetores() antes de compute_JV(0)\n";
    compute_JV(0);
      cout<< "PhElem<NumVariaveis>::inicia_vetores() apos de compute_JV(0)\n";
  }
};

// ************************************
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::finaliza_vetores()
{
  if(vetores_iniciados == 1) {
    vetores_iniciados = 0;
    //cout<< "PhElem<NumVariaveis>::finaliza_vetores()\n";
    delete [] JV; JV=nullptr;
    //delete [] Jb; Jb=nullptr;

    int k;
    for (k=0;k<NumVariaveis;k++){
      delete [] u0[k]; u0[k]=nullptr;
      //   delete [] ua[k];ua[k]=nullptr;
      delete [] usave[k];usave[k]=nullptr;
    }
    delete [] stgbtrbmapM;stgbtrbmapM=nullptr;
    delete [] stgbtrbmapP;stgbtrbmapP=nullptr;
  }
};

#endif /* PhElem_hpp */
