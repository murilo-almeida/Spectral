#ifndef _DGProb_headers
#define _DGProb_headers
// ****************************************************************************
// Classe DG : Derivada de GeProb (Generic Problem)
// Encapsula os dados para resolucao do problema algebrico generico
// ****************************************************************************

#include "PhElem.hpp"
#include "GeProb.hpp"
#include "DG_Elem.h"
#include "Fluids.h"

// ************************************************************
// Esta linha eh importante: MyElem eh usado em varios lugares
typedef DG_Elem  MyElem;
// observe que: class DG_Elem : public PhElem<2>
// ************************************************************

// GeProb com elemento MyElem=PhElem<2> = e duas variaveis;
// o segundo 2 abaixo equivale a dois espacos interpolantes

class DG_Prob : public GeProb<MyElem,2,2>
{
 public:

  DG_Prob(Epetra_Comm& comm);
  ~DG_Prob();
  void preamble(char * str);
  void Driver(char * argv = 0);

  void DG_alocar_mem_local(const int qmax,const int nsat, const int npres);
  void DG_liberar_mem_local(const int nsat, const int npres);
  void DG_initial_conditions();

  // void DG_eco(std::vector< std::vector <int> > MapRow, int * ncount);
  void DG_eco();
  void DG_Iterate();
  void DG_Iterate_NOX();
  void Ge_MVRA0(const double Dt,
                Epetra_Map Map,
                double & valor0,
                double_t & norm_delta_X);
  void Ge_MVRA(const double Dt,
               Epetra_Map Map,
               double & valor,
               int & token,
               double_t & norm_delta_X);

  int Eigenvectors(const double Dt, const Epetra_Map & Map);
  void DG_conditionNumber(Epetra_Map Map);
  void Ler_Arquivo_Dados_DG(char *arq_par);
  void MPI_Recebe_Dados_DG();
  void DG_Transfere_rst(int * b_in,double * b_do,double * buff);
  void DG_MatrizVetor_UMFPACK(const double Dt,int & count, const int nz,
                      int * Ti, int * Tj, double * Tx, double * B);
  void DG_MatrizVetor_Epetra(const double Dt,
                             Teuchos::RCP<Epetra_FECrsMatrix> A,
                             Teuchos::RCP<Epetra_FEVector> RHS);
  void DG_EI_Interior(const EDGE border, int & count,
                      int * Ti,int * Tj, double * Tx,
                      double * B);
  void DG_EI_Inflow(const EDGE border, int & count,
                    int * Ti,int * Tj, double * Tx,
                    double * B);
  void DG_EI_Outflow(const EDGE border,  int & count,
                     int * Ti,int * Tj, double * Tx,
                     double * B);
  void DG_EI_Epetra_Interior(const EDGE border,
                             Teuchos::RCP<Epetra_FECrsMatrix> A,
                             Teuchos::RCP<Epetra_FEVector> RHS);
  void DG_EI_Epetra_Inflow(const EDGE border,
                           Teuchos::RCP<Epetra_FECrsMatrix> A,
                           Teuchos::RCP<Epetra_FEVector> RHS);
  void DG_EI_Epetra_Outflow(const EDGE border,
                      Teuchos::RCP<Epetra_FECrsMatrix> A,
                      Teuchos::RCP<Epetra_FEVector> RHS);
  int DG_EI_Inflow(const EDGE border,
                   double * mx,
                   double * B,
                   int    * indx);
  int DG_EI_Interior(const EDGE border,
                     double * mx,
                     double * B,
                     int    * indx);
  int DG_EI_Outflow(const EDGE border,
                    double * mx,
                    double * B,
                    int    * indx);
  /*void DG_EI(Epetra_FECrsMatrix & A, Epetra_FEVector & RHS); */
  void DG_EI_flux(const EDGE border, double & flux_w, double & flux_n);
  // Initial Guess Edge Integrals
  void DG_EI_IG(const EDGE border,
                const double sigma,
                const double sigma1,
                const double beta,int & count,
                int * Ti,int * Tj, double * Tx,
                double * B,
                const int qmax,
                const double w0[]);
  // ****************************
  //void RowEI(const EDGE border, int * NumNz, int ** Map);
  //void ValoresInterpolados(char * str = 0);
  void DG_Escrever_rst(const int nprt=0);
	void DG_imprimir_taxas_de_producao(FILE * fout3,
                                           double valor, double valor0, double valor1,
                                           int iter);

  void projetar_C0(FILE *file,double (*func)(double,double,double),
                   const int & ivar);

  bool evaluate(const FillType flag,
                const Epetra_Vector* X,
                Epetra_Vector* FVec,
                Epetra_RowMatrix* Jacobian);
  
  void Processa_condicoes_contorno(/*int *BC, std::vector< std::vector<int> > face_mask*/);
  // ***************************************************************************

 private:

  // const int NumVARS = 2  //!< especificados em GeProb<2,2>
  // const int NumFIELDS = 2 //!< especificados em GeProb<2,2>
  //int P0;//!< ordem dos polinômios da variável 0 (sat)
  //int P1;//!< ordem dos polinômios da variável 1 (pres)
  //int  Q;//!< numero de pontos de Gauss em cada direcao
  // int nsat;
  //int npres;
  int count;
  int nz;
  int restart_flag;//!< 0 comeca; 1 continua
  double theta,pd,mun,muw,rhon,rhow;
  double porosidade;
  char arq_rst[256];//!< string contendo o nome do arquivo para recomecar simulacao
  char arq_sai[256];//!< string contendo o nome do arquivo para saida dos resultados
  char arq_eco[256];//!< string contendo o nome do arquivo para eco dos dados
  char arq_flu[256];//!< string contendo o nome do arquivo de fluxos
  double t;//!< tempo de simulacao
  double Dt;//!< Intervalo de tempo; pode variar conforme falhe a convergencia
  double Dt0;//!< Intervalo de tempo inicial
  double fatordt; //!< fator de reducao do Dt
  double t_final;//!< tempo final de simulacao
  double tprint; //!< intervalo entre escritas de resultados no arquivo
  double precisao;//!< precisao da iteracao de Gauss
  double qdado;
  int passo;
  int tnc;
  int tncut;
  //double p_in, p_out, sn_in, sn_ini, pw_ini;
	double injrate_w,injrate_n;
  double prodrate_w,prodrate_n;
  double Iw,In,Qw,Qn;
  double injrate_w0,injrate_n0;
  double prodrate_w0,prodrate_n0;
  double gravidade[3];
  double perm ;
  double sigma,sigma1 ;
  double beta  ;

  Fluids fluids;

  // **********************
  // alocar memoria local *
  // **********************
  // Tracos de pw e sn; acompanham os edges
  double * gbtrpw;
  double * gbtrsn;
  //double * phi_r;
  //double * phi_l;
  double *** Kgsn;
  double *** Kgpw;
  double *** Kgpc;
  /*!< \short Kgxx [iE][i][iq]
   * = Traco de Kgrad_xx
   * \param  iE  = Elemento da borda = 0 ou 1 (lado da borda)
   * \param  dir = direcao da derivada
   * \param  iq = indice do ponto de quadratura de Gauss
   */

 // double ** gphi_r;
 // double ** gphi_l;
  //double ** gphi_ ;
  /*!< \short gphi_# [dir][iq]
   * = grad_phi
   * \param  dir = direcao da derivada
   * \param  iq = indice do ponto de quadratura de Gauss
   */

  // Armazenar trphi para [elemento][variavel][modo][iq]
  // double **** trphi;
  // /*!< \short trphi [elemento][variavel][modo][iq]
  //  * = traco de phi
  // * \param Elemento da borda = 0 ou 1 (lado da borda)
  // * \param Variavel = sat (=0) ou pres (=1)
  // * \param  Modo = numero do modo local
  // * \param  iq = indice do ponto de quadratura de Gauss
  // */

  double **** K_g_phi_n;
  /*!< \short K_g_phi_n [iE][variavel][modo][iq]
   * = Traco K * grad_phi . normal
   * \param iE = Elemento da borda = 0 ou 1 (lado da borda)
   * \param Variavel = sat (=0) ou pres (=1)
   * \param  Modo = numero do modo local
   * \param  iq = indice do ponto de quadratura de Gauss
   */
  bool memoria_alocada = false;
};

#endif
