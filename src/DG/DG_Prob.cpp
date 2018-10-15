
/*****************************************************************************/
/*****************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DG_run"
/*****************************************************************************/
/*****************************************************************************/
#include "DG_Prob.h"
#include "PhElem.hpp"
#include "GeProb.hpp"

// Definicao da classe DG_Prob
// Deriva de GeProb
DG_Prob::DG_Prob(Epetra_Comm& comm) : GeProb<MyElem,2,2>::GeProb(comm)
{};

DG_Prob::~DG_Prob()
{
	cout << "\nDG_Prob::~DG_Prob()\n";
	// *************************************************************************
	// Liberar memoria local
	// *************************************************************************
	const int nsat=el[0].show_ptr_stdel(sat)->nn_val();// Num de modos para saturacao
	const int npres=el[0].show_ptr_stdel(pres)->nn_val();// Num de modos para pressao
	cout << "Ponto1 nsat ="<< nsat<< " npres = "<< npres << "\n";
	DG_liberar_mem_local(nsat,npres);

	for(int i = 0; i < NELEM; i++) {
		el[i].finaliza_vetores();
	}


		delete [] gbtrpw;gbtrpw=nullptr;
		delete [] gbtrsn;gbtrsn=nullptr;

		cout << "Final de DG_Prob::~DG_Prob()\n";
	//*/
};


/*****************************************************************************/
/*****************************************************************************/
/*
void  DG_Prob::set_orders()
{
  int n = Field.size();
  ptrLinear.set(n);
  ptrTriang.set(n);
  ptrQuadri.set(n);

  for(int i=0; i < n; ++i) {
  //  cout << "DG_Prob::set_orders:  ";
  // cout << "Aloca memoria dinamica\n ";
  ptrLinear[i] = new        Linear(Field[i].P,Field[i].Q);
  ptrTriang[i] = new      Triangle(Field[i].P,Field[i].Q);
  ptrQuadri[i] = new Quadrilateral(Field[i].P,Field[i].Q);
  }
};
 */
// ****************************************************************************
// ****************************************************************************
void DG_Prob::DG_initial_conditions()
{
  FILE * frestart;
  // printf("comeco de DG_Prob::DG_initial_conditions\n");

  fluids.set_mun(mun);
  fluids.set_muw(muw);
  fluids.set_rhon(rhon);
  fluids.set_rhow(rhow);
  fluids.set_theta(theta);
  fluids.set_pd(pd);

  NumD=0;

	//FILE * f_eco;
	// f_eco = fopen("saida_eco","wb");

	for(int i = 0; i < NELEM; i++) {
    //printf("ponto 0 :Elemento %d\n",i);
    el[i].set_porosidade(porosidade);
    el[i].set_permeabilidade(perm,perm,perm);
    el[i].set_fontes(0.0,0.0);

    // *********************************************************************
    // mapeia os modos locais nos globais; calcula NumD;
    // Eh importante; nao pode ser descartado
    // NumD eh usado para definir o tamanho dos vetores B, x e xc abaixo
    // *********************************************************************
    el[i].inicia_gbnmap(NumD);

    //Inicaliza os vetores locais: JV,b0,bs,Grad's e tracos
    el[i].inicia_vetores();

		// elemento necessita ter pontos de Gauss nas bordas
    //el[i].inicia_tracos(border); // calcula tracos usando pontos de Gauss nas bordas

		// elemento nao necessita ter pontos de Gauss nas bordas
		el[i].inicia_funcoes_na_borda(border); // para calcular tracos diretamente nas bordas

		//fprintf(f_eco,"elemento %d\n",i);
    //el[i].echo_traco(f_eco);
  }

	//fclose(f_eco);

  // ***************************************************************************
  // Exemplo
  // caso com fontes qq9.neu
  // fprintf(fout1,"\nFontes de sw no elemento 0\n\n");
  //el[0].set_fontes(qdado,0.0);// qw=qdado; qn=0.0;
  //el[80].set_fontes(0.0,-qdado);
  // *****************************************************************************

  printf("Numero de variaveis globais NumD: %d\n",NumD);


  const int ndim = el[0].show_ptr_stdel(0)->ndim_val();// dimensao espacial do dominio
  // const int NGQP = el[0].show_ptr_stdel(0)->NGQP_val();// numero de pontos de Gauss
  const int qmax=el[0].show_ptr_stdel(0)->qborder_val();// Pontos de gauss nas arestas
  const int nsat=el[0].show_ptr_stdel(sat)->nn_val();// Num de modos para saturacao
  const int npres=el[0].show_ptr_stdel(pres)->nn_val();// Num de modos para pressao

  // alocar memoria local para arrays de DG_Prob para conter tracos
  DG_alocar_mem_local(qmax,nsat,npres);

  if(restart_flag==0) {// inicio da simulacao

    int nz=10000000; // Ponto critico; nz pequeno pode dar segmentation fault
    Tstruct T(nz);
		Vstruct B(NumD);
		Vstruct x(NumD);

    t=0.0;

    // *************************************
    // Montar as condicoes iniciais        *
    // *************************************
    // inicializa a variavel sn

    for(int i=0;i<NELEM;i++) {
      el[i].projetar_C0(NULL,sn_inicial,sat);
    }

    //fprintf(fout1,"\n");
    // for(int i=0;i<NELEM;i++) {
    // //caso com fontes e sem fluxo nas laterais
    // el[i].projetar_C0(fout1,pw_inicial,pres);
    //  }
    // fprintf(fout1,"\n\n");

    // inicializa a variavel sn
    NumD = 0;// contador de incognitas de sn

    //  gbnmap_discontinuous(sat,NumD);
    // Loop abaixo eh equivalente ao gbnmap_discontinuous() acima
    for(int i=0;i<NELEM;i++)
      el[i].inicia_gbnmap(sat,NumD);
    printf("\nNumero de incognitas de saturacao: NumD = %d\n",NumD);


    // refaz gbnmap para resolver soh pw
    NumD = 0;// contador de incognitas de pw

		//gbnmap_continuous(pres,NumD); // pressao continua ao longo das bordas // 13/12/2013

    // pressao descontinua ao longo das bordas
    // gbnmap_discontinuous(pres,NumD);
    for(int i=0;i<NELEM;i++)
      el[i].inicia_gbnmap(pres,NumD);
    printf("Numero de incognitas de pressao: NumD = %d\n",NumD);


    // Montar a matriz e vetor

    int count=0;
    // inicializa o vetor B
    for(int i=0;i<NumD;i++)
      B.v[i]=0.0;
    // Integral de volume
    for(int i=0;i<NELEM;i++) {
      //  cout << "VolumeIntegrals_IG elemento "<< i << endl;
      el[i].VolumeIntegrals_IG(fluids,count,T.i,T.j,T.x,B.v);
    }
    printf("integral de volume count= %d\n",count);
    // Integral sobre os borders
    double sigma0=1.0e10;

		printf ("qmax = %d\n",qmax);

    DG_EI_IG(sigma0,sigma0,beta,count,T.i,T.j,T.x,B.v);

    if(nz<count) exit(0);
    //fprintf(fout1,"apos integrais: count = %d  nz = %d\n",count,nz);
    printf("myid = %d apos integrais: count = %d  nz = %d\n",myid,count,nz);

    /*
    //Trecho para imprimir vetor e matriz para usar em outro programa
    FILE * fmatriz;
    fmatriz=fopen("matriz.dat","wb");
    fprintf(fmatriz,"%d %d\n",NumD,count);
    for(int i=0;i<NumD;i++)fprintf(fmatriz,"%g\n",B[i]);
    for(int i=0;i<count;i++)
    fprintf(fmatriz,"%d %d %g\n",Ti[i],Tj[i],Tx[i]);
    fclose(fmatriz); //exit(0);
    */

    //resolve o problema DG_IG para achar o pw0
    // Usa UMFPACK
    ResolverSistema(NumD,count,T.i,T.j,T.x,B.v,x.v);
    printf(" Resolveu sistema\n");

    // repassar os resultados para os elementos
    for(int i=0;i<NELEM;i++) {
      // cout << " Elemento "<< i << endl;
      el[i].mapa_inverso(pres,x.v); // Copia os coeficientes de x em u0
    }

    // atualizar os valores de sna e pwa e escrever as variaveis em fout
    if(myid==0) {
      FILE *fout;
      fout=fopen(arq_sai,"wb");
      fprintf(fout,"# t = %e\n",t);
      for(int i=0;i<NELEM;i++)
        el[i].Atualizar_valores(fout);// atualiza os valores das concentracoes no pontos de Gauss
      fprintf(fout,"\n\n");
      fclose(fout);
    }
    else {
      for(int i=0;i<NELEM;i++) {
        el[i].Atualizar_valores();
        // atualiza os valores das concentracoes no pontos de Gauss
      }
    }
    printf("Atualizou valores\n");
    // *******************************************
    // Preparar condicao inicial de fluxos
    // *******************************************
    // Valores de vazao
    Iw = 0.0;
    In = 0.0;
    Qn = 0.0;
    Qw = 0.0;
    injrate_w0 = 0.0;
    injrate_n0 = 0.0;
    prodrate_w0 = 0.0;
    prodrate_n0 = 0.0;
    double auxw, auxn;
    // printf("inj_w0 %e inj_n0 %e\n",injrate_w0,injrate_n0);
    //printf("pro_w0 %e pro_n0 %e\n",prodrate_w0,prodrate_n0);
    for(int ii = 0; ii < nin; ii++) {
      DG_EI_flux(border[in_borders[ii]],auxw,auxn);
      injrate_w0 += auxw;
      injrate_n0 += auxn;
    }
    //printf("inj_w0 %e inj_n0 %e\n",auxw,auxn);
    for(int ii = 0; ii < nout; ii++) {
      DG_EI_flux(border[out_borders[ii]],auxw,auxn);
      prodrate_w0 -= auxw;
      prodrate_n0 -= auxn;
    }
    // printf("myid = %d prod_w0 %e prod_n0 %e\n",myid,prodrate_w0,prodrate_n0);

		Dt=Dt0;

  }// fim de restart_flag==0
  else { // continuar : restart_flag==1
    // *************************************************************************
    // Continuar simulacao
    // *************************************************************************
    // Ler  frestart
		double  buff [NELEM*nsat*npres];
    int     b_in [3];
    double  b_do [10];
    int count1=0;
    if(myid==0) {
      frestart=fopen(arq_rst,"rb");
      fscanf(frestart,"%d %d",&b_in[0],&b_in[1]);
      for(int nn=0;nn<10;nn++)
        fscanf(frestart,"%lf",&b_do[nn]);
      for(int i=0;i<NELEM;i++)
        el[i].ler_restart_buffer(frestart,buff,count1);
      // for(int i=0;i<NELEM;i++)el[i].ler_restart(frestart);
      fclose(frestart);
      b_in[2]=count1;
    }

#ifdef HAVE_MPI
    //MPI::COMM_WORLD.Barrier();
    Comm->Barrier();
    DG_Transfere_rst(b_in,b_do,buff);
#endif

    // atualizar os valores de sna e pwa
    // cout << "Atualizar valores inicial \n"<< endl;
    for(int i=0;i<NELEM;i++)
      el[i].Atualizar_valores();// atualiza os valores das concentracoes no pontos de Gauss
    // cout << "SAIU Atualizar valores inicial \n"<< endl;
  }

  // ***********************************************
  // refazer gbnmap para incluir as duas variaveis
  // ***********************************************

	 // pressao descontinua ao longo das bordas // 13/12/2013
	  NumD=0; // necessita ser especificado para renumerar todas as variaveis
#ifdef GBNMAP_DG
  printf("Usa GBNMAP_DG : Elementos descontinuos\n");
  // fprintf(fout1,"Usa _DG_GBNMAP : Elementos descontinuos\n");
  for(int i=0;i<NELEM;i++) {
    el[i].inicia_gbnmap(NumD);// mapeia os modos locais nos globais
  }                           // para as duas variaveis
#else
  printf("Usa gbnmap_continuous_2D : Elementos continuos\n");
  //fprintf(fout1,"Usa gbnmap_continuous_2D : Elementos continuos\n");

  gbnmap_continuous(NumD); // Todas as variaveis sao continuas Usada em 22/10/2016
#endif

  // Mantem a pressao continua ao longo das bordas // 13/12/2013
	//for(int i=0;i<NELEM;i++) {
  //  el[i].inicia_gbnmap(sat,NumD);// mapeia os modos locais nos globais
	// }
	// 13/12/2013

  // *************************************************************************
  //cout << "Preparou configuracao inicial\n";
  cout << "Final de DG_Prob::DG_initial_conditions " << "myid =" << myid << "\n";
}// Fim de DG_initial_conditions()


// *****************************************************************************
// *****************************************************************************
// Constroi a Matriz e o Vetor para usar com UMFPACK
// Esta funcao e chamada so por DG_eco para dimensionar a matriz e vetor
// *****************************************************************************
void DG_Prob::DG_MatrizVetor_UMFPACK(const double Dt,int & count, const int nz,
                                     int * Ti, int * Tj, double * Tx, double * B)
{
  cout << "DG_MatrizVetor_UMFPACK\n ";
  EDGE border1;
  //  montar matriz e vetor
  // inicializar o vetor B
  for(int i=0;i<NumD;i++)
    B[i]=0.0;
  count=0;

  for (int p=0;p<NumPart;p++) {
    // elementos da particao
    for (int k=0;k<Particao[p].nele;k++) {
      int i=Particao[p].ele[k];
      el[i].VolumeIntegrals_UMFPACK(Dt,fluids,count,Ti,Tj,Tx,B,gbtrsn,gbtrpw);
      //printf("terminou VolumeIntegralsUMFPACK %2d  count= %2d\n",i,count);
    }
    // Elementos ghosts
    for (int k=0;k<Particao[p].ngho;k++) {
      int i=Particao[p].gho[k];
      el[i].VolumeTracos(Dt,fluids,gbtrsn,gbtrpw);
    }

    for (int k=0;k<Particao[p].nbor;k++) {
      // bordas da particao
      int i=Particao[p].bor[k];

      border1=border[i];
      int t=border1.tipo;

      switch(t) {

      case -1:
					DG_EI_Inflow(border1,count,Ti,Tj,Tx,B);
					break;

				case 1:
					DG_EI_Outflow(border1,count,Ti,Tj,Tx,B);
					break;

				case 2:
					DG_EI_Interior(border1,count,Ti,Tj,Tx,B);
					break;
      }
    }
  }

  //printf("count na saida do switch %d \n",count);
  if(nz<count) {
    printf("nz<count= %3d\n",count);
    exit(0);
  }
};

// ****************************************************************************
// ****************************************************************************
void DG_Prob::DG_alocar_mem_local(const int qmax,
                                  const int nsat,
                                  const int npres)
{

  // **********************
  // alocar memoria local *
  // **********************

  // Para armazenar tracos;
  // Nao confundir com variaveis locais em VolumeIntegrals_IG e DG_VI
	// Que sao da classe PhElem
  //phi_r = new double [qmax];
  //phi_l = new double [qmax];

  Kgsn = new double ** [2];
  Kgpw = new double ** [2];
  Kgpc = new double ** [2];
  for (int i=0; i<2; i++) {
    Kgsn[i]= new double * [2];
    Kgpw[i]= new double * [2];
    Kgpc[i]= new double * [2];

    for(int j=0; j<2; j++) {
      Kgsn[i][j]= new double [qmax];
      Kgpw[i][j]= new double [qmax];
      Kgpc[i][j]= new double [qmax];
    }
  }

  K_g_phi_n = new double *** [2];
  for(int iE=0;iE<2;iE++) {
    K_g_phi_n[iE]= new double ** [2];

    K_g_phi_n[iE][sat]=new double * [nsat];
    K_g_phi_n[iE][pres]=new double * [npres];
    for(int iM=0;iM<nsat;iM++) {
      K_g_phi_n[iE][sat][iM] = new double [qmax];
    }
    for(int iM=0;iM<npres;iM++) {
      K_g_phi_n[iE][pres][iM] = new double [qmax];
    }
  }
  memoria_alocada = true;
  cout << "DG_Prob::DG_alocar_mem_local : Sucesso!\n";
};

// *****************************************************************
// *****************************************************************
void DG_Prob::DG_liberar_mem_local(const int nsat, const int npres)
{
  if(memoria_alocada) {
	// *************************
	// liberar a memoria local *
	// *************************

//	delete []  in_borders;  in_borders=nullptr;
	//delete [] out_borders; out_borders=nullptr;

	for (int i=0; i<2; i++) {
		for(int j=0; j<2; j++) {
			delete [] Kgsn[i][j]; Kgsn[i][j]=nullptr;
			delete [] Kgpw[i][j]; Kgpw[i][j]=nullptr;
			delete [] Kgpc[i][j]; Kgpc[i][j]=nullptr;
		}
		delete [] Kgsn[i];Kgsn[i]=nullptr;
		delete [] Kgpw[i];Kgpw[i]=nullptr;
		delete [] Kgpc[i];Kgpc[i]=nullptr;
		//delete [] gphi_r[i];gphi_r[i]=nullptr;
		//delete [] gphi_l[i];gphi_l[i]=nullptr;
		//delete [] gphi_[i]; gphi_[i] =nullptr;
	}
	delete [] Kgsn;Kgsn=nullptr;
	delete [] Kgpw;Kgpw=nullptr;
	delete [] Kgpc;Kgpc=nullptr;
	//delete [] gphi_r; gphi_r=nullptr;
	//delete [] gphi_l; gphi_l=nullptr;
	//delete [] gphi_ ; gphi_ =nullptr;

	//delete [] phi_r;phi_r=nullptr;
	//delete [] phi_l;phi_l=nullptr;

	// Liberar memoria de trphi e K_g_phi_n
	for(int iE=0;iE<2;iE++) {
		for(int iM=0;iM<nsat;iM++) {
			//  delete [] trphi[iE][sat][iM];trphi[iE][sat][iM]=0;
			delete []  K_g_phi_n[iE][sat][iM]; K_g_phi_n[iE][sat][iM]=nullptr;
		}
		for(int iM=0;iM<npres;iM++) {
			//delete [] trphi[iE][pres][iM]; trphi[iE][pres][iM]=0;
			delete []  K_g_phi_n[iE][pres][iM]; K_g_phi_n[iE][pres][iM]=nullptr;
		}
		//delete [] trphi[iE][sat]; trphi[iE][sat]=0;
		//delete [] trphi[iE][pres];trphi[iE][pres]=0;
		//delete [] trphi[iE];trphi[iE]=0;
		delete [] K_g_phi_n[iE][sat];K_g_phi_n[iE][sat]=nullptr;
		delete [] K_g_phi_n[iE][pres]; K_g_phi_n[iE][pres]=nullptr;
		delete [] K_g_phi_n[iE]; K_g_phi_n[iE]=nullptr;
	}
	//delete [] trphi; trphi=nullptr;
	delete [] K_g_phi_n;K_g_phi_n=nullptr;
  memoria_alocada = false;
	cout << "DG_Prob::DG_liberar_mem_local : Sucesso!\n";
}
};

// ************************************************************************
// Escrever arquivo restart
// ************************************************************************
void DG_Prob::DG_Escrever_rst(const int nprt)
{
  FILE * frestart;
  char string[256];
  sprintf(string,"%s.%d",arq_sai,nprt);
  frestart=fopen(string,"wb");
  fprintf(frestart,"%d %d %e %e\n",passo,tnc,Dt,t);
  fprintf(frestart,"%e %e %e %e ",injrate_w0,injrate_n0,Iw,In);
  fprintf(frestart,"%e %e %e %e\n",prodrate_w0,prodrate_n0,Qw,Qn);
  for(int i=0;i<NELEM;i++)
    el[i].escrever_restart(frestart);
  fclose(frestart);
};
// ************************************************************************
void DG_Prob::projetar_C0(FILE *file,double (*func)(double,double,double),
												 const int & ivar)
{
  //FILE * frestart;
  printf("\n\ncomeco de DG_Prob::projetar_C0\n");
  /*
  fluids.set_mun(mun);
  fluids.set_muw(muw);
  fluids.set_rhon(rhon);
  fluids.set_rhow(rhow);
  fluids.set_theta(theta);
  fluids.set_pd(pd);
  */
  NumD=0;

	//FILE * f_eco;
	// f_eco = fopen("saida_eco","wb");

	for(int i = 0; i < NELEM; i++) {
    //printf("ponto 0 :Elemento %d\n",i);
    el[i].set_porosidade(porosidade);
    el[i].set_permeabilidade(perm,perm,perm);
    el[i].set_fontes(0.0,0.0);

    // *********************************************************************
    // mapeia os modos locais nos globais; calcula NumD;
    // Eh importante; nao pode ser descartado
    // NumD eh usado para definir o tamanho dos vetores B, x e xc abaixo
    // *********************************************************************
    el[i].inicia_gbnmap(NumD);

    //Inicaliza os vetores locais: JV,b0,bs,Grad's e tracos
    el[i].inicia_vetores();

		// elemento necessita ter pontos de Gauss nas bordas
    //el[i].inicia_tracos(border); // calcula tracos usando pontos de Gauss nas bordas

		// elemento nao necessita ter pontos de Gauss nas bordas
		//el[i].inicia_funcoes_na_borda(border); // para calcular tracos diretamente nas bordas

		//fprintf(f_eco,"elemento %d\n",i);
    //el[i].echo_traco(f_eco);
  }

	//fclose(f_eco);

  // ***************************************************************************
  // Exemplo
  // caso com fontes qq9.neu
  // fprintf(fout1,"\nFontes de sw no elemento 0\n\n");
  //el[0].set_fontes(qdado,0.0);// qw=qdado; qn=0.0;
  //el[80].set_fontes(0.0,-qdado);
  // *****************************************************************************

  printf("DG_Prob::projetar_C0: Numero de variaveis globais NumD: %d\n",NumD);

 // const int ndim = el[0].show_ptr_stdel(0)->ndim_val();// dimensao espacial do dominio
  // const int NGQP = el[0].show_ptr_stdel(0)->NGQP_val();// numero de pontos de Gauss
  const int qmax=el[0].show_ptr_stdel(0)->qborder_val();// Pontos de gauss nas arestas
  const int nsat=el[0].show_ptr_stdel(sat)->nn_val();// Num de modos para saturacao
  const int npres=el[0].show_ptr_stdel(pres)->nn_val();// Num de modos para pressao

  // alocar memoria local para arrays de DG_Prob para conter tracos
  DG_alocar_mem_local(qmax,nsat,npres);
  /*
  if(restart_flag==0) {// inicio da simulacao

    int nz=5000000; // Ponto critico; nz pequeno pode dar segmentation fault
    Tstruct T(nz);
		Vstruct B(NumD);
		Vstruct x(NumD);

    t=0.0;

    // *************************************
    // Montar as condicoes iniciais        *
    // *************************************
    // inicializa a variavel sn
   */
    for(int i=0;i<NELEM;i++) {
     // cout << "\nChamar PhElem::projetar_C0 para elemento "<< i << endl;
      el[i].projetar_C0(file,func,ivar);
     // cout << "Terminou PhElem::projetar_C0 para elemento "<< i << endl<< endl;
    }

};
// ******************************************************************************************
