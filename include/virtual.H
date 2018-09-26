virtual void set(int P,int Q);
virtual void Jacobian(const Vertice vert[],const int map[], double * JV );

virtual double mass(int, int, const double []);
virtual void printtofile(FILE * fout,const double u[], 
			 double (*)(double,double,double),
			 const Vertice vert[], const int map[]);
virtual void printtofile(FILE * fout,const double u[],  
			 const Vertice vert[], const int map[]);
virtual void printGQtofile(FILE * fout,const double ftemp[],
			   const double ftemp1[],
			   const Vertice vert[], const int map[]);
virtual void printwGQtofile(FILE * fout,
			    const Vertice vert[],
			    const int map[],
			    const double JV[]) ;
virtual void computeVertice(double f_vert[],const double u[], 
			    const Vertice vert[], const int map[]);
virtual void computeAtPoints(const int npoints, const double LocCoord[],const double u[],
			     const Vertice vert[], const int map[],double f[],double GloCoord[]);
virtual void evalGQ(double f0[],double f1[],
		    const double u0[],const double u1[]);
virtual void evalGQ(double f0[],const double u0[],
                    const int NF=1,
                    const int nvar=0);
virtual void eval_Phi(const int m,double Phi[]);
virtual void eval_GradPhi(const Vertice vert[], const int map[],const int m,double ** der);
/*
virtual void Processar_geometria(int nel,const Vertice * vert,
				 const int numv,
				 const int * VN,
				 int map[],
				 int sgn[],int sinal_temp[],
				 int& NG,int& NL,
				 std::vector<EDGE>& border,
				 int Ng[],
				 int & NF, int Face[], int Fng[],
				 int f_mask[],
				 const int N_add=0);
*/
virtual void Gradiente(FILE * fout, double * grad[],
		       double (*func)(double, double, double),
		       const Vertice vert[], const int map[]);
virtual void Gradiente(FILE * fout, double * grad[],
		       const double fvec[],
		       const Vertice vert[], const int map[]);
virtual void Gradiente(double * grad[],
		       const double fvec[],
		       const Vertice vert[], const int map[]);
virtual void print_nome(FILE *);

virtual void Dirichlet(const int aresta,
		       const Vertice vert[],const int vert_map[],
		       const int nmap[],const int sgn[],
		       int bflag[],double Xbc[],
		       double (*f)(double,double,double));
//   virtual void BoundForce(const int aresta,const Vertice vert[],
// 			  const int gbnmap[],double bf[],
// 			  double (*f)(double,double,double),
// 			  const int varn,const int NFields);
//   
//   virtual void fluxo(const int aresta,const Vertice vert[],
// 		     const int gbnmap[],const double f0[][3],
// 		     double bs[],
// 		     const int varn,const int NFields);
virtual void teste(int& );
virtual void trace(const int lado,const int Q,const int sinal,
		   const double *valores,double *saida);
virtual void computeFuncGQ(double f_[], 
			   const Vertice vert[], const int map[],
			   double (*func)(double,double,double));
virtual void vector_of_integral_of_f_Phi_dv(double vec[],
				  const double func[],
				 // const Vertice vert[],const int map[],
				  const double JV[]);
virtual void vector_of_integral_of_f_Phi_dv(double vec[],
				  double (*func)(double, double, double),
				  const Vertice vert[], const int map[],
				  const double JV[]);
virtual void elem_traces(const Vertice ptvert[],const int Vert_map[],const int sinal[],
			 double *** TP,double **** TGP, double * Jb);
virtual void trace_Jb(const Vertice vert[],const int map[],const int sinal[],
                      double * Jb);
virtual const int aresta_lvert(const int & i, const int & j) const;
virtual const int face_lvert(const int & i, const int & j) const;
virtual const int show_nvf(const int &i) const;
virtual const int show_face_tipo(const int &i) const;
virtual const int show_fd0(const int &i) const;
virtual const int show_fd1(const int &i) const;
virtual const int show_fv2(const int &i) const;
virtual const int show_ind_mode(const int & i, const int & j = 0, const int & k =0) const;
virtual void superficie_externa(const int Vert_map[], const Vertice vert[],
                                const int & num_local,double & area,double normal[3]);

private:

//virtual void make_mass(double M[], const double JV[]);
virtual void make_mass_matrices(int nfields);
virtual void printtoarray(const double u[], 
			  const Vertice vert[], const int map[],
			  double x[], double y[], double z[],
			  double ftemp[]);
//  virtual void Ler_geometria(int nel,const Vertice * vert,FILE * finput,
// 			     int map[],int sgn[],
// 			     int& NG,int& NL,EDGE border[],int Ng[], 
// 			     int & NF, int Face[], int Fng[], int f_mask[]);

//virtual void set_b(Mat2<double> & b, const Vertice vert[], const int map[]){};
