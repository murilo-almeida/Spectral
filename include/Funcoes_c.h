#ifndef _Funcoes_c_headers
#define _Funcoes_c_headers

// ****************************************************************************
void jacobi_poly ( int n, double alpha, double beta, double x, double cx[] );
// ****************************************************************************
// The 3 polymorphic forms of Jacobi_P
// ****************************************************************************
void Jacobi_P (int n, double a, double b, double x, 
	       double &y, double &dy, double &d2y);
void Jacobi_P (int n, double a, double b, double x, 
	       double &y, double &dy);
void Jacobi_P (int n, double a, double b, double x, 
	       double &y);
// ****************************************************************************

// ****************************************************************************
// Raizes dos polinomios de Jacobi de grau m (alpha, beta)
// ****************************************************************************
void Jacobi_roots(int m, double alpha, double beta, double x[]);
// ****************************************************************************

// ****************************************************************************
// Calculo da funcao gamma
// ****************************************************************************
double gamma_log ( double x );
double my_gamma ( double x );
// ****************************************************************************

// ****************************************************************************
// Funcoes envolvidas no calculo dos pontos e pesos de integracao
// de Gauss_Jacobi
// ****************************************************************************
void Gauss_Jacobi_parameters(int Q, double alpha, double beta, 
				     double x[], double w[],double D[MAXQ][MAXQ]);

// ****************************************************************************
// Funcoes envolvidas no calculo dos pontos e pesos de integracao
// de Gauss_Lobatto_Jacobi
// ****************************************************************************
void Gauss_Lobatto_Jacobi_parameters(int Q, double alpha, double beta, 
				     double x[], double w[],double D[MAXQ][MAXQ]);
double C(int Q2, double alpha, double beta, double xi);
// *******************************************************************

// ****************************************************************************
// Funcoes envolvidas no calculo dos pontos e pesos de integracao             *
// de Gauss_Radau_Jacobi                                                      *
// ****************************************************************************
void Gauss_Radau_Jacobi_parameters(int Q, double alpha, double beta,
				   double x[], double w[],double D[MAXQ][MAXQ]);
double B(int q, double alpha, double beta, double xi);
// ****************************************************************************

double d_epsilon ( void );

// ****************************************************************************
// Funcoes Principais Modificadas                                             *
// ****************************************************************************
double Psia(int I, int i, double z);
double Psib(int I, int J, int i, int j, double z);
double Psic(int I, int J, int K, int i, int j, int k, double z);
// ****************************************************************************
// Derivadas                                                                  *
// ****************************************************************************
double DPsia(int I, int i, double z);
double DPsib(int I, int J, int i, int j, double z);
double DPsic(int I, int J, int K, int i, int j, int k, double z);
// ****************************************************************************
// Collocation Differentiation                                                *
// Gauss-Lobato-Jacobi, Gauss-Radau-Jacobi e Gauss-Jacobi                     *
// ****************************************************************************
void ColDifGLJ(int Q,double alpha,double beta,double D[MAXQ][MAXQ]);
void ColDifGRJ(int Q,double alpha,double beta,double D[MAXQ][MAXQ]);
void ColDifGJ (int Q,double alpha,double beta,double D[MAXQ][MAXQ]);


// ****************************************************************************
// Funcoes globais para especifidar as condicoes externas                     *
// ****************************************************************************
//void teste_aresta(int p,int n0,int n1,int& sgn, int& ng,
	//	  int& NG,int& NL,int Na[],int Nb[],int Ng[]);
int teste_aresta(int p,int n0,int n1,int& sgn, int& ng,
                 int& NG,int& NL, std::vector<EDGE>& border,int Ng[],int nel,int nlado,
		  const Vertice * vert,const int sinal_normal);
int aresta_gbnum(const int na, const int nb, int & NL, std::vector<ARESTA> & aresta_vector);
int face_gbnum(const int nvf, const int v[], int & NF, std::vector<FACE> & face_vector);
double funcao_contorno(double, double,double);
double funcao(double, double, double);
double g(double, double, double);
double force0(double, double, double);
double force1(double, double, double);
void predictor(int N,double *fn3,double *fn2,double *fn1,double *fn,double*yn,
	       double h, double *yp);

//void imprime(FILE * fout, double x, double * y,
//	     int NELEM,int NG,int NumD, int * novoNum, PhElem el[],
//	     double X0[], double X1[]);

int get_linha(FILE * finput,char linha[256]);


double p1(double,double,double);
double p2(double,double,double);
double fluxo_entrando(double,double,double);
double funcao_pdir(double,double,const int);
double funcao_sdir(double,double);
double sn_inicial(double x, double y,double z);
double pw_inicial(double x, double y,double z);
double norma_max(int N,double V[]);
double norma_l2(int N,double B[]);
void ordenar(const int n,int base[],int seq[]);
void quad_ordem(const int Av[4],const int Bi[4],int dir[3], int sgn[3], int & v2);
void tetrahedro_faz_face_mask(const int & n_in,int ver_temp[],std::vector<int> & face_mask);

// ****************************************************************************
void somavv(int q,const double * ptx,const double * pty,double * ptr);
void soma3v(int q,const double *x1,const double *x2,const double *x3,
	    double * ptr);
void soma4v(int q,const double *x1,const double *x2,const double *x3,
	    const double *x4,double * ptr);
void soma5v(int q,const double *x1,const double *x2,const double *x3,
	    const double *x4,const double *x5,
	    double * ptr);
void subtvv(int q,const double * ptx,const double * pty,double * ptr);
void prodvv(int q,const double * ptx,const double * pty,double * ptr);
void prod3v(int q,const double * ptx,const double * pty, const double * ptz,
	    double * ptr);
void prod4v(int q,const double *x1,const double *x2, const double *x3,
	    const double *x4,double * ptr);
void prod5v(int q,const double *x1,const double *x2,const double *x3,
	    const double *x4,const double *x5,
	    double * ptr);
void prodev(int q, const double fator, const double * ptx, double *ptr);
void prodgg(int q,const double *ptx0, const double *ptx1,const double *pty0,
	    const double *pty1, 
	    double *ptr);
void prodgn(int q,const double *ptx0, const double *ptx1, const double *n,
	    double *ptr);
void prodgn(const int ndim,const int q,double **grad, const double *n,
	    double *ptr);
void Kgradn(int q,const double *K,const double *ptx0, const double *ptx1, 
	    const double *n,double *ptr);
void Kgradn(const int ndim,const int qmax,const double *K,double **grad, 
	    const double *n,double *ptr);
void Kgradgrad(int q,const double *K,const double *ptx0, const double *ptx1, 
	       const double *pty0, 
	       const double *pty1,double *ptr);
void integral(int q,const double *ptw, const double *ptx, double  &r);
void soma_vetor(int q,const double *ptx, double & r);
// ****************************************************************************
int Sinal(int i);
int ResolverSistema(const int NumD, const int count,
		    int * Ti, int * Tj, double * Tx,
		    double * B, double *x);
int CH_ResolverSistema(const int NumD, const int count,
		       int * Ti, int * Tj, double * Tx,
		       double * B, double *x);


//void teste_tracoII(FILE * fout1,
//									 const Vertice vert[],
//									 const PhElem el[],
//									 const EDGE border);

#endif
