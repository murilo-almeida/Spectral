#ifndef _Stdel_headers
#define _Stdel_headers
// ****************************************************************************
class Stdel
// ****************************************************************************
{ 
 public:
  
  Stdel() {};
  ~Stdel() {};
  int show_vtk_type() const {return vtk_type;};
  void printStdel() const;
  int tipo_val() const {return tipo;};
  int ndim_val() const {return ndim;};
  int qborder_val() const {return qborder;};
  int nborder_val() const {return nborder;};
  int NGQP_val() const {return NGQP;};
  int nv_val() const {return nv;};
  int ne_val() const {return ne;};
  int nf_val() const {return nf;};
  int nn_val() const {return nn;};
  int nb_val() const {return nb;};
  int P_val(int i) const {return P[i];};
  int Q_val(int i) const {return Q[i];};
  int gqt_val(int i) const {return gqt[i];};
  //double w_val(int i,int j) const {return wGQ[i][j];};
 // Mode ind_val(int i) const {return mode_[i];};
  double show_Phi_val(const int i, const int j)const {return Phi_val[i][j];};

  double show_Mb_comp(int, int) const; 
  double show_Mi_inv(int,int) const;
  double show_McMi_inv(int,int) const;
  void print_matrices(FILE *);
  void show_ind(int i,int& p,int& q,int& r);
  void duplicar_mass_matrix(int Nfields);
  int w_vec(double saida[]);
  int show_emapv(int i)const { return emapv[i];};
  int show_emapi(int i)const { return emapi[i];};
  int show_bmapv(int i)const { return bmapv[i];};
  int show_bmapi(int i)const { return bmapi[i];};
  
  const int is_on_border(int m, int e,int & pos) const
  {
      int aux=0;
      int ans=0; // NO
      for(int i=emapi[e]; (i < emapi[e+1] && ans == 0); i++)
        {
            if(emapv[i] == m){
                ans=1; // YES
                pos=aux;
            }
            aux++;
        }
      return ans;
  };
  
  void print_border_map(FILE * fout)
  {
    for(int i=0;i<ne; i++){
      fprintf(fout,"Edge %d\n",i);
      for(int j=emapi[i];j<emapi[i+1];j++)
	fprintf(fout,"%d ",emapv[j]);
      fprintf(fout,"\n");
    }
  };
  
  // Pure virtual functions
  virtual void set(int P,int Q)=0;
  virtual void Jacobian(const Vertice vert[],const int map[], double * JV )=0;
  
  virtual double mass(int, int, const double [])=0;
  virtual void printtofile(FILE * fout,const double u[], 
                           double (*)(double,double,double),
                           const Vertice vert[], const int map[])=0;
  virtual void printtofile(FILE * fout,const double u[],  
                           const Vertice vert[], const int map[])=0;
  virtual void printGQtofile(FILE * fout,const double ftemp[],
                             const double ftemp1[],
                             const Vertice vert[], const int map[])=0;
  virtual void printwGQtofile(FILE * fout,
                              const Vertice vert[],
                              const int map[],
                              const double JV[]) =0;
  virtual void computeVertice(double f_vert[],const double u[], 
			      const Vertice vert[], const int map[])=0;
  virtual void computeAtPoints(const int npoints, const double LocCoord[],const double u[],
			       const Vertice vert[], const int map[],double f[],double GloCoord[])=0;
  virtual void evalGQ(double f0[],double f1[],const double u0[],const double u1[])=0;
  virtual void evalGQ(double f0[],const double u0[],const int NF=1,const int nvar=0)=0;
  virtual void eval_Phi(const int m,double Phi[])=0;
  virtual void eval_GradPhi(const Vertice vert[], const int map[],const int m,double ** der)=0;
  virtual void Gradiente(FILE * fout, double ** grad,  
			 double (*func)(double, double, double), 
			 const Vertice vert[], const int map[])=0;
  virtual void Gradiente(FILE * fout, double ** grad,  
			 const double fvec[], 
			 const Vertice vert[], const int map[])=0;
  virtual void Gradiente(double ** grad,  
			 const double fvec[], 
			 const Vertice vert[], const int map[])=0;
  virtual void print_nome(FILE *)=0;

  virtual void Dirichlet(const int aresta,
			 const Vertice vert[],const int vert_map[],
			 const int mmap[],const int sgn[],
			 int bflag[],double Xbc[],
			 double (*f)(double,double,double))=0;
    virtual void face_Jacobian(const int face_num,
                               const Vertice vert[],
                               const int vert_map[], // numero global dos vertices dos nos
                               const int sgn[],
                               double * J)=0;
//   virtual void BoundForce(const int aresta,const Vertice vert[],
// 			  const int gbnmap[],double bf[],
// 			  double (*f)(double,double,double),
// 			  const int varn,const int NFields)=0;
//   
//   virtual void fluxo(const int aresta,const Vertice vert[],
// 		     const int gbnmap[],const double f0[][3],
// 		     double bs[],
// 		     const int varn,const int NFields)=0;
  virtual void teste(int& )=0;
  // Calcular os tracos nos pontos de quadratura Gauss_Jacobi (nao incluir os pontos extremos
  virtual void trace(const int lado,const int Q,const int sinal,
                     const double *valores,double *saida, const int vert_map[] = 0)=0;
  virtual void computeFuncGQ(double f_[], 
			     const Vertice vert[], const int map[],
			     double (*func)(double,double,double))=0;
  virtual void vector_of_integral_of_f_Phi_dv(double vec[],
				    const double func[],
				    //const Vertice vert[],const int map[],
				    const double JV[])=0;
  virtual void vector_of_integral_of_f_Phi_dv(double vec[],
				    double (*func)(double, double, double), 
				    const Vertice vert[], const int map[],
				    const double JV[])=0; 
  virtual void elem_traces(const Vertice ptvert[],const int Vert_map[],
                           const int sinal[],double *** TP,double **** TGP, double * Jb)=0;
  virtual void trace_Jb(const Vertice vert[],const int map[],const int sinal[],
                        double * Jb)=0;
  virtual const int aresta_lvert(const int & i, const int & j) const = 0;
  virtual const int face_lvert(const int & i, const int & j) const = 0;
  virtual const int show_nvf(const int &i) const = 0;
  virtual const int show_face_tipo(const int &i) const = 0;
  virtual const int show_fd0(const int &i) const = 0;
  virtual const int show_fd1(const int &i) const = 0;
  virtual const int show_fv2(const int &i) const = 0;
  virtual const int show_ind_mode(const int & i, const int & j = 0, const int & k =0)
                                  const = 0;
  virtual void superficie_externa(const Vertice vert[],const int Vert_map[], 
                                  const int & num_local,double & area,double normal[3]) = 0;
  // virtual void set_b(Mat2<double> & b, const Vertice vert[], const int map[])=0;
 // virtual const int aresta_val(int i, int j) = 0;// const {return aresta[i][j];};

  
 protected:  
  
  int tipo;    ///< type of the element
  int ndim;    ///< number of geometrical dimension 
  int nv;      ///< number of vertices
  int ne;      ///< number of edges
  int nf;      ///< number of faces
  int nborder; ///< number of borders
  int qborder;    ///< number of Gauss quadrature points on the faces (borders)
  int vtk_type;///< number of element for vtk file
  int nn;      ///< number of modes (and/or modes)
  int nb;      ///< number of boundary (shared) modes (or modes)
  int P[3];    ///< (number of modes - 1) in each direction 
  int Q[3];    ///< number of Gauss Quadrature points in direction i
  int NGQP;    ///< total number of Gauss Quadrature points
  int gqt[3];  /*!< gqt Gauss Quadrature Type:
                  \param      1= Gauss-Jacobi
                  \param      2= Gauss-Radau-Jacobi
                  \param      3= Gauss-Lobatto-Jacobi
                  \param Index of the row refers to the direction:
                  \param 0=x;
                  \param 1=y;
                  \param 2=z.
                  */
                
  double xGQ[3][MAXQ];/*!< coordinates of the Gauss quadrature points
		\param Index of the row refers to the direction: 0=x; 1=y; 2=z.
                \param Index of the column refers to the quadrature point */
  double wGQ[3][MAXQ];  /*!< weights of the Gauss quadrature points
	\param Index of the row refers to the direction: 0=x; 1=y; 2=z.
                \param Index of the column refers to the quadrature point */
  Mode mode_[MAXMODES];///< mode_[m(ijk)]= number of the expansion function into the local
  
  double D[MAXQ][MAXQ][3]; /*!< Coeficientes para a derivacao por colocacao 
			     Matrizes locais */
  double Phi_val[MAXMODES][MAXQ*MAXQ]; ///< Contem os vetores Phi nos pontos de quadratura
  double Mb_comp[MAXNBL][MAXNBL]; //!< Matriz de massa compacta de Schur
  double Mi_inv[MAXNIL][MAXNIL];//!< Matriz de massa inversa dos modos internos de Schur
  double McMi_inv[MAXNBL][MAXNIL];//!< Matriz de massa inversa dos modos compacta de Schur
  int * emapv= nullptr; //!< edge_map_valores
  int * emapi= nullptr; //!< edge_map_inicio
  int * bmapv= nullptr; //!< border_map_valores
  int * bmapi= nullptr; //!< border_map_inicio
  double *  D_Phi_val = nullptr;
  
  //virtual void make_mass_matrices(int NFIELDS)=0;
// virtual void printtoarray(const double u[], 
// 			    const Vertice vert[], const int map[],
// 			    double x[], double y[], double z[], 
// 			    double ftemp[])=0;

};

#endif
